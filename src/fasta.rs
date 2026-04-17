/// A parsed FASTA record.
pub struct Record {
    /// Header (first word after `>`, no `>` prefix)
    pub header: String,
    /// Sequence data, uppercased, concatenated from all lines
    pub seq: Vec<u8>,
}

impl Record {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Parse FASTA from a reader. Mirrors the behavior of MaxBin2's fastaReader::parse():
/// - Header is truncated at first whitespace
/// - Sequence is uppercased
/// - Empty sequences are included (matching original behavior)
pub fn parse<R: std::io::BufRead>(reader: R) -> Vec<Record> {
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = Vec::new();

    // Matches fastaReader.cpp:523-581 (parse() main loop)
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => break,
        };
        // Matches fastaReader.cpp:527-530: strip trailing \n and \r
        let line = line.trim_end_matches(['\n', '\r']).to_string();

        if let Some(rest) = line.strip_prefix('>') {
            // Matches fastaReader.cpp:531-556: header line handling
            if let Some(header) = current_header.take() {
                // Matches fastaReader.cpp:538-539: addSeq() called for previous record
                records.push(Record {
                    header,
                    seq: current_seq,
                });
                current_seq = Vec::new();
            }
            // Matches fastaReader.cpp:548-555: truncate header at first space or tab
            let header = rest.split_whitespace().next().unwrap_or("").to_string();
            current_header = Some(header);
        } else if current_header.is_some() && !line.is_empty() {
            // Matches fastaReader.cpp:562-571: accumulate sequence lines
            // Matches fastaReader.cpp:488-499 (convertToUpper): lowercase to uppercase
            current_seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
        // Matches fastaReader.cpp:558-560: lines before first '>' are silently skipped
    }

    // Matches fastaReader.cpp:582-586: flush last record after EOF
    if let Some(header) = current_header {
        records.push(Record {
            header,
            seq: current_seq,
        });
    }

    records
}

/// Parse FASTA from a file path. Handles .gz transparently.
pub fn parse_file(path: &std::path::Path) -> std::io::Result<Vec<Record>> {
    let file = std::fs::File::open(path)?;

    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        let decoder = flate2::read::GzDecoder::new(file);
        let reader = std::io::BufReader::new(decoder);
        Ok(parse(reader))
    } else {
        let reader = std::io::BufReader::new(file);
        Ok(parse(reader))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn basic_parse() {
        let input = b">seq1 some description\nACGT\nTGCA\n>seq2\naaaa\n";
        let records = parse(Cursor::new(input));

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "seq1");
        assert_eq!(records[0].seq, b"ACGTTGCA");
        assert_eq!(records[0].len(), 8);
        assert_eq!(records[1].header, "seq2");
        assert_eq!(records[1].seq, b"AAAA"); // uppercased
    }

    #[test]
    fn header_truncated_at_whitespace() {
        let input = b">contig_1\tflag=1 multi=4.98 len=2926\nACGT\n";
        let records = parse(Cursor::new(input));
        assert_eq!(records[0].header, "contig_1");
    }

    #[test]
    fn empty_sequence() {
        let input = b">empty\n>notempty\nACGT\n";
        let records = parse(Cursor::new(input));
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "empty");
        assert_eq!(records[0].seq, b"");
        assert_eq!(records[1].seq, b"ACGT");
    }

    #[test]
    fn skips_lines_before_first_header() {
        let input = b"some junk\nmore junk\n>real\nACGT\n";
        let records = parse(Cursor::new(input));
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, "real");
    }

    #[test]
    fn handles_crlf() {
        let input = b">seq1\r\nACGT\r\nTGCA\r\n";
        let records = parse(Cursor::new(input));
        assert_eq!(records[0].seq, b"ACGTTGCA");
    }

    /// Decompress a .gz file to a temp path for the C++ parser (which can't read gzip).
    fn decompress_gz(gz_path: &std::path::Path) -> std::path::PathBuf {
        use std::io::{Read, Write};
        let out_path = std::env::temp_dir().join("maxbin_rs_test_contigs.fa");
        let file = std::fs::File::open(gz_path).unwrap();
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut contents = Vec::new();
        decoder.read_to_end(&mut contents).unwrap();
        let mut out = std::fs::File::create(&out_path).unwrap();
        out.write_all(&contents).unwrap();
        out_path
    }

    #[test]
    fn equivalence_with_original_fasta_parser() {
        // This test requires nf-core test data (only available in devshell).
        let gz_path = match std::env::var("MAXBIN2_TEST_CONTIGS") {
            Ok(p) => std::path::PathBuf::from(p),
            Err(_) => return,
        };

        // Parse with our Rust implementation (handles gzip natively)
        let rust_records = parse_file(&gz_path).expect("rust parser failed");

        // Decompress for the C++ parser, then parse with the original
        let plain_path = decompress_gz(&gz_path);
        let original = crate::original_ffi::OriginalFastaReader::new(&plain_path);

        // Same number of records
        assert_eq!(
            rust_records.len(),
            original.num_records() as usize,
            "record count mismatch"
        );

        // Field-by-field comparison for every record
        for (i, rust_rec) in rust_records.iter().enumerate() {
            let idx = i as u32;

            assert_eq!(
                rust_rec.header,
                original.header(idx),
                "header mismatch at record {i}"
            );

            assert_eq!(
                rust_rec.seq,
                original.seq(idx),
                "sequence mismatch at record {i} (header: {})",
                rust_rec.header
            );

            assert_eq!(
                rust_rec.len(),
                original.seq_len(idx) as usize,
                "length mismatch at record {i} (header: {})",
                rust_rec.header
            );
        }

        std::fs::remove_file(&plain_path).ok();
    }

    #[test]
    fn parse_nfcore_test_data() {
        // This test uses the nf-core test contigs if available via env var.
        // Verifies we parse the same number of records as the original MaxBin2.
        let path = match std::env::var("MAXBIN2_TEST_CONTIGS") {
            Ok(p) => std::path::PathBuf::from(p),
            Err(_) => return, // skip outside devshell
        };

        let records = parse_file(&path).expect("failed to parse test contigs");

        // Verified against: zcat test1.contigs.fa.gz | grep -c "^>" => 272
        assert_eq!(records.len(), 272, "expected 272 contigs");

        for r in &records {
            assert!(!r.header.is_empty(), "empty header found");
        }

        // MaxBin2 default min_contig_length=1000. From the reference run:
        // 219 contigs pass, 53 go to .tooshort
        let long_contigs = records.iter().filter(|r| r.len() >= 1000).count();
        assert_eq!(long_contigs, 219);
        assert_eq!(records.len() - long_contigs, 53);
    }
}
