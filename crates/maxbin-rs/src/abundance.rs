/// A parsed abundance record: contig header -> abundance value.
pub struct AbundanceRecord {
    pub header: String,
    pub abundance: f64,
}

/// Minimum abundance value — matches VERY_SMALL_NUM in the original.
/// Matches AbundanceLoader.cpp:3 (`#define VERY_SMALL_NUM 0.0001`)
const VERY_SMALL_NUM: f64 = 0.0001;

/// Parse an abundance file. Mirrors MaxBin2's AbundanceLoader::parse() behavior:
/// - Lines are split at the first tab, space, comma, or semicolon
/// - Header has trailing whitespace stripped
/// - Abundance values below VERY_SMALL_NUM are clamped to VERY_SMALL_NUM
/// - Empty lines are skipped
/// - BUG (reproduced): the "skip multiple separators" loop on line 115 of
///   AbundanceLoader.cpp uses && instead of ||, so it never executes.
///   We reproduce this: only the first separator character is consumed.
pub fn parse<R: std::io::BufRead>(reader: R) -> Result<Vec<AbundanceRecord>, String> {
    let mut records = Vec::new();

    // Matches AbundanceLoader.cpp:84-139 (parse() main loop)
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => return Err(format!("read error: {e}")),
        };

        // Matches AbundanceLoader.cpp:90-93: strip trailing \n, \r, space, tab
        let line = line.trim_end_matches(['\n', '\r', ' ', '\t']);

        // Matches AbundanceLoader.cpp:94-97: skip empty lines
        if line.is_empty() {
            continue;
        }

        // Matches run_MaxBin.pl:1521-1523 (checkAbundFile): strip leading '>'
        // from headers. Some abundance files use FASTA-style ">contig_name"
        // headers; the Perl preprocesses these before the C++ sees them.
        let line = line.strip_prefix('>').unwrap_or(line);

        // Matches AbundanceLoader.cpp:99-107: advance past the header until separator
        let sep_pos = line
            .find(['\t', ' ', ',', ';'])
            .ok_or_else(|| format!("no separator found in line: {line}"))?;

        let header_part = &line[..sep_pos];

        // The original has a bug here (AbundanceLoader.cpp:115):
        // `while (*c == '\t' && *c == ' ' && *c == ',' && *c == ';')` uses
        // && instead of ||, so multiple consecutive separators are never
        // skipped. We fix this by trimming all leading separators/whitespace.
        let value_part = line[sep_pos..].trim_start_matches(['\t', ' ', ',', ';']);

        if value_part.is_empty() {
            return Err(format!("no value after separator in line: {line}"));
        }

        // Matches AbundanceLoader.cpp:164-168: trim trailing spaces/tabs from header
        let header = header_part.trim_end().to_string();

        let abundance: f64 = value_part.parse().unwrap_or(0.0);

        // Matches AbundanceLoader.cpp:173-176: clamp to VERY_SMALL_NUM
        let abundance = if abundance < VERY_SMALL_NUM {
            VERY_SMALL_NUM
        } else {
            abundance
        };

        records.push(AbundanceRecord { header, abundance });
    }

    Ok(records)
}

/// Parse an abundance file from a path.
pub fn parse_file(path: &std::path::Path) -> Result<Vec<AbundanceRecord>, String> {
    let file =
        std::fs::File::open(path).map_err(|e| format!("can't open {}: {e}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    parse(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn basic_parse() {
        let input = b"contig1\t42.5\ncontig2\t0.001\n";
        let records = parse(Cursor::new(input)).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "contig1");
        assert_eq!(records[0].abundance, 42.5);
        assert_eq!(records[1].header, "contig2");
        assert_eq!(records[1].abundance, 0.001);
    }

    #[test]
    fn clamps_very_small_values() {
        let input = b"contig1\t0.00001\ncontig2\t0.0\n";
        let records = parse(Cursor::new(input)).unwrap();

        assert_eq!(records[0].abundance, VERY_SMALL_NUM);
        assert_eq!(records[1].abundance, VERY_SMALL_NUM);
    }

    #[test]
    fn skips_empty_lines() {
        let input = b"\ncontig1\t1.0\n\ncontig2\t2.0\n\n";
        let records = parse(Cursor::new(input)).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn accepts_various_separators() {
        let tab = parse(Cursor::new(b"c1\t1.0\n" as &[u8])).unwrap();
        let space = parse(Cursor::new(b"c1 1.0\n" as &[u8])).unwrap();
        let comma = parse(Cursor::new(b"c1,1.0\n" as &[u8])).unwrap();
        let semi = parse(Cursor::new(b"c1;1.0\n" as &[u8])).unwrap();

        assert_eq!(tab[0].abundance, 1.0);
        assert_eq!(space[0].abundance, 1.0);
        assert_eq!(comma[0].abundance, 1.0);
        assert_eq!(semi[0].abundance, 1.0);
    }

    #[test]
    fn strips_fasta_header_prefix() {
        let input = b">contig1\t42.5\n>contig2\t0.001\ncontig3\t1.0\n";
        let records = parse(Cursor::new(input)).unwrap();

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].header, "contig1");
        assert_eq!(records[1].header, "contig2");
        assert_eq!(records[2].header, "contig3");
    }

    #[test]
    fn multiple_separators_skipped() {
        // The original has a bug where multiple consecutive separators are
        // not skipped (&&  instead of ||). We fix this — multiple tabs,
        // spaces, etc. between header and value are handled correctly.
        let input = b"c1\t\t1.0\nc2\t  \t2.5\n";
        let records = parse(Cursor::new(input)).unwrap();

        assert_eq!(records[0].header, "c1");
        assert_eq!(records[0].abundance, 1.0);
        assert_eq!(records[1].header, "c2");
        assert_eq!(records[1].abundance, 2.5);
    }

    #[test]
    fn equivalence_with_original() {
        // Writes a temp abundance file, parses with both Rust and C++ FFI,
        // compares field by field.
        use std::io::Write;

        let test_data = "contig_1\t10.5\ncontig_2\t0.00001\ncontig_3\t100.0\ncontig_4 0.5\n";

        let tmp_path = std::env::temp_dir().join("maxbin_rs_test_abund.txt");
        {
            let mut f = std::fs::File::create(&tmp_path).unwrap();
            f.write_all(test_data.as_bytes()).unwrap();
        }

        // Rust
        let rust_records = parse_file(&tmp_path).unwrap();

        // Original C++
        let original = crate::original_ffi::OriginalAbundanceLoader::new(&tmp_path);

        assert_eq!(rust_records.len(), original.num_records() as usize);
        assert!(original.is_parse_success());

        for (i, rec) in rust_records.iter().enumerate() {
            let orig_abund = original.abundance_by_index(i as i32);
            assert!(
                (rec.abundance - orig_abund).abs() < 1e-10,
                "abundance mismatch at record {i} (header: {}): rust={} original={}",
                rec.header,
                rec.abundance,
                orig_abund,
            );
        }

        std::fs::remove_file(&tmp_path).ok();
    }
}
