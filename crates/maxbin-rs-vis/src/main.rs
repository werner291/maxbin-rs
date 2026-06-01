//! maxbin-rs visualizer.
//!
//! Section 1 - "What are the inputs". The EM stage does not see raw reads or
//! the assembly; it sees the *pipeline intermediaries* left by the earlier
//! stages: a contigs FASTA, one or more abundance (depth) profiles, and a seed
//! list of initial bin centers. Pick a dataset, commit it, and the page reports
//! what those intermediaries actually contain. Nothing happens until you
//! commit a choice.
//!
//! Everything stays in the page: the intermediaries are fetched as static
//! assets and parsed client-side. No server, no upload.

mod copy;
mod data;

use data::Inputs;
use leptos::prelude::*;
use leptos::task::spawn_local;

fn main() {
    console_error_panic_hook::set_once();
    leptos::mount::mount_to_body(App);
}

/// The dataset choices offered in the dropdown. `key` selects the loader arm.
const DATASETS: &[(&str, &str)] = &[
    ("bfragilis", "B. fragilis (single genome, one sample, smoke test)"),
    ("cami-i-high", "CAMI I High (5-sample benchmark metagenome)"),
];

#[component]
fn App() -> impl IntoView {
    // The dropdown selection, not yet committed.
    let choice = RwSignal::new(String::new());
    // The committed dataset's parsed metadata; None until "Select to continue".
    let inputs = RwSignal::new(None::<Inputs>);
    let status = RwSignal::new(String::new());
    let loading = RwSignal::new(false);

    let commit = move |_| {
        let key = choice.get();
        if key.is_empty() || loading.get() {
            return;
        }
        loading.set(true);
        inputs.set(None);
        status.set(format!("Loading {key} intermediaries…"));
        spawn_local(async move {
            match load(&key).await {
                Ok(i) => {
                    status.set(String::new());
                    inputs.set(Some(i));
                }
                Err(e) => status.set(format!("Failed: {e}")),
            }
            loading.set(false);
        });
    };

    view! {
        <div class="field">
            <main class="sheet">
                <p class="kicker">{copy::KICKER}</p>
                <h1>{copy::PAGE_TITLE}</h1>
                <p class="intro">{rich(copy::INTRO_1)}</p>
                <p class="intro">{rich(copy::INTRO_2)}</p>

                <h2>{copy::TITLE}</h2>
                <p class="lede">{copy::LEDE_1}</p>
                <p class="intro">{copy::LEDE_2}</p>
                <p class="intro">{copy::LEDE_3}</p>

                <div class="chooser">
                    <select
                        prop:value=move || choice.get()
                        on:change=move |ev| choice.set(event_target_value(&ev))
                    >
                        <option value="">"choose a dataset…"</option>
                        {DATASETS.iter().map(|(k, label)| view! {
                            <option value=*k>{*label}</option>
                        }).collect_view()}
                    </select>
                    <button
                        on:click=commit
                        disabled=move || choice.get().is_empty() || loading.get()
                    >
                        {move || if loading.get() { "Loading…" } else { "Select to continue →" }}
                    </button>
                </div>

                {move || (!status.get().is_empty()).then(|| view! {
                    <p class="status">{status.get()}</p>
                })}

                {move || inputs.get().map(metadata_view)}
            </main>
        </div>
    }
}

/// Render the metadata for a committed dataset.
fn metadata_view(i: Inputs) -> impl IntoView {
    let contigs = i.contigs.clone();
    let samples = i.samples.clone();
    view! {
        <section class="meta">
            <h2>{i.label}</h2>
            <p class="blurb">{i.blurb}</p>

            <h3>{copy::PIPELINE_HEADING}</h3>
            <p class="note">{copy::PIPELINE_INTRO}</p>
            {pipeline_view(&i)}

            <h3>"Contigs"</h3>
            {match contigs {
                Some(c) => {
                    let gc = c.gc.map(|g| format!("{:.1}%", g * 100.0))
                        .unwrap_or_else(|| "n/a".into());
                    view! {
                        <table>
                            <tbody>
                                <tr><th>"count"</th><td class="num">{group(c.n as u64)}</td></tr>
                                <tr><th>"total length"</th><td class="num">{format!("{} bp", group(c.total_bp))}</td></tr>
                                <tr><th>"length min / median / max"</th>
                                    <td class="num">{format!("{} / {} / {} bp", group(c.len_min as u64), group(c.len_med as u64), group(c.len_max as u64))}</td></tr>
                                <tr><th>"GC content"</th><td class="num">{gc}</td></tr>
                            </tbody>
                        </table>
                    }.into_any()
                }
                None => view! { <p class="note">{i.contigs_note.clone()}</p> }.into_any(),
            }}

            <h3>{format!("Abundance: {} ({} contigs profiled)", i.abund_format, group(i.abund_rows as u64))}</h3>
            <table>
                <thead><tr><th>"sample"</th><th class="num">"mean depth"</th><th class="num">"min"</th><th class="num">"max"</th></tr></thead>
                <tbody>
                    {samples.into_iter().map(|s| view! {
                        <tr>
                            <th>{s.label}</th>
                            <td class="num">{format!("{:.2}×", s.mean)}</td>
                            <td class="num">{format!("{:.2}", s.min)}</td>
                            <td class="num">{format!("{:.2}", s.max)}</td>
                        </tr>
                    }).collect_view()}
                </tbody>
            </table>

            <h3>"Seeds"</h3>
            {match i.seed_n {
                Some(n) => view! { <p>{format!("{n} seed contigs: {}.", copy::SEED_DESC)}</p> }.into_any(),
                None => view! { <p class="note">{copy::NO_SEED_NOTE}</p> }.into_any(),
            }}

            {(!i.notes.is_empty()).then(|| view! {
                <div class="notes">
                    {i.notes.into_iter().map(|n| view! { <p class="note">{n}</p> }).collect_view()}
                </div>
            })}
        </section>
    }
}

/// The pipeline chain, as a vertical rail, with each step's per-dataset status.
/// The EM is the only step that runs in the browser; the rest are external
/// tools whose outputs each dataset may or may not bundle.
fn pipeline_view(i: &Inputs) -> impl IntoView + use<> {
    let has_abundance = !i.samples.is_empty();
    let has_seeds = i.seed_n.is_some();
    let runnable = i.has_sequence && has_abundance && has_seeds;

    // What the EM is still missing, if anything.
    let mut missing: Vec<&str> = Vec::new();
    if !i.has_sequence {
        missing.push("contig sequence");
    }
    if !has_abundance {
        missing.push("abundance");
    }
    if !has_seeds {
        missing.push("seeds");
    }

    // (css-class, name, tool · io, chip-class, chip-text)
    let em_chip = if runnable {
        ("live", "runs live in your browser".to_string())
    } else {
        ("missing", format!("blocked: needs {}", missing.join(" + ")))
    };
    let abund_chip = if has_abundance {
        ("have", "abundance provided")
    } else {
        ("missing", "abundance not bundled")
    };
    let seed_chip = if has_seeds {
        ("have", "seeds provided")
    } else {
        ("missing", "seeds not bundled")
    };

    let rows: Vec<(&str, &str, &str, &str, String)> = vec![
        (
            "upstream",
            "Assemble",
            "assembler (not part of maxbin) · reads → contigs",
            "upstream",
            "happens before maxbin".to_string(),
        ),
        (
            "offline",
            "Map reads, count depth",
            "bowtie2 · contigs + reads → abundance",
            abund_chip.0,
            abund_chip.1.to_string(),
        ),
        (
            "offline",
            "Find single-copy marker genes",
            "FragGeneScanRs + HMMER · contigs → seeds",
            seed_chip.0,
            seed_chip.1.to_string(),
        ),
        (
            "browser",
            "Bin by expectation-maximization",
            "maxbin-rs · tetranucleotides + abundance + seeds → bins",
            em_chip.0,
            em_chip.1,
        ),
    ];

    view! {
        <div class="pipeline">
            {rows.into_iter().map(|(run, name, tool, chip, chip_text)| view! {
                <div class=format!("stage {run}")>
                    <div class="stage-head">
                        <span class="stage-name">{name}</span>
                        <span class=format!("chip {chip}")>{chip_text}</span>
                    </div>
                    <div class="stage-tool">{tool}</div>
                </div>
            }).collect_view()}
        </div>
    }
}

/// Render a string with `*word*` spans as emphasis; everything else plain.
/// Odd-indexed segments (between asterisks) are the emphasized ones.
fn rich(text: &'static str) -> impl IntoView {
    text.split('*')
        .enumerate()
        .map(|(i, seg)| {
            if i % 2 == 1 {
                view! { <em>{seg}</em> }.into_any()
            } else {
                view! { {seg} }.into_any()
            }
        })
        .collect_view()
}

/// Thousands-separated integer (e.g. 2_803_653_876 -> "2,803,653,876").
fn group(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    let bytes = s.as_bytes();
    for (idx, b) in bytes.iter().enumerate() {
        if idx > 0 && (bytes.len() - idx) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

/// Fetch and parse the intermediaries for one dataset.
async fn load(key: &str) -> Result<Inputs, String> {
    match key {
        "bfragilis" => {
            let base = "assets/bfragilis";
            let contigs = get(&format!("{base}/contigs.fa")).await?;
            let abund = get(&format!("{base}/abund")).await?;
            let seed = get(&format!("{base}/seed")).await.unwrap_or_default();

            let stats = data::fasta_stats(&contigs);
            let (rows, sample) = data::maxbin_abund(&abund);
            Ok(Inputs {
                label: "B. fragilis",
                blurb: copy::BFRAGILIS_BLURB,
                contigs: Some(stats),
                has_sequence: true,
                contigs_note: String::new(),
                abund_format: "MaxBin 2-column",
                abund_rows: rows,
                samples: vec![sample],
                seed_n: Some(data::seed_count(&seed)),
                notes: vec![],
            })
        }
        "cami-i-high" => {
            // The MetaBAT depth file is the real abundance intermediary and
            // carries each contig's length, so we get contig stats from it
            // without shipping the 2.8 GB assembly FASTA.
            let depth = get("assets/cami-i-high/depth").await?;
            let (stats, samples, rows) = data::metabat_depth(&depth);
            Ok(Inputs {
                label: "CAMI I High-Complexity",
                blurb: copy::CAMI_BLURB,
                contigs: Some(stats),
                has_sequence: false,
                contigs_note: String::new(),
                abund_format: "MetaBAT depth",
                abund_rows: rows,
                samples,
                seed_n: None,
                notes: vec![copy::CAMI_NOTE.to_string()],
            })
        }
        other => Err(format!("unknown dataset '{other}'")),
    }
}

/// GET a static asset as text.
async fn get(url: &str) -> Result<String, String> {
    gloo_net::http::Request::get(url)
        .send()
        .await
        .map_err(|e| e.to_string())?
        .text()
        .await
        .map_err(|e| e.to_string())
}
