//! maxbin-rs visualizer — EM-input stage.
//!
//! Loads the EM subcommand's three inputs (contigs FASTA, abundance, seed),
//! parses them client-side, and scatter-plots the contigs (abundance vs length,
//! seeds highlighted). No EM run yet — that lands once `maxbin-core` is
//! wasm-extracted. Everything stays in the page: built-in datasets are fetched
//! as static assets, custom uploads are read with the File API, no server.

mod data;

use data::{Contig, Dataset};
use leptos::prelude::*;
use leptos::task::spawn_local;
use std::f64::consts::TAU;
use wasm_bindgen::JsCast;
use web_sys::{CanvasRenderingContext2d, HtmlCanvasElement, HtmlInputElement};

fn main() {
    console_error_panic_hook::set_once();
    leptos::mount::mount_to_body(App);
}

const CANVAS_W: u32 = 640;
const CANVAS_H: u32 = 420;

#[component]
fn App() -> impl IntoView {
    let dataset = RwSignal::new(None::<Dataset>);
    let status = RwSignal::new("Pick a dataset.".to_string());
    let custom = RwSignal::new(false);
    let canvas_ref = NodeRef::<leptos::html::Canvas>::new();

    // Load a built-in dataset by fetching its three asset files in parallel-ish.
    let load_builtin = move |name: &'static str| {
        status.set(format!("Loading {name}…"));
        spawn_local(async move {
            let base = format!("assets/{name}");
            match fetch_triple(&base).await {
                Ok((c, a, s)) => {
                    let ds = Dataset::from_inputs(&c, &a, &s);
                    status.set(format!("{}: {} contigs", name, ds.contigs.len()));
                    dataset.set(Some(ds));
                }
                Err(e) => status.set(format!("Failed to load {name}: {e}")),
            }
        });
    };

    // Repaint whenever the dataset changes (or the canvas mounts).
    Effect::new(move |_| {
        let ds = dataset.get();
        let Some(canvas) = canvas_ref.get() else {
            return;
        };
        let canvas: HtmlCanvasElement = canvas;
        let ctx = canvas
            .get_context("2d")
            .unwrap()
            .unwrap()
            .dyn_into::<CanvasRenderingContext2d>()
            .unwrap();
        draw(&ctx, ds.as_ref());
    });

    // Custom upload: three file inputs read into these signals, then joined.
    let contigs_text = RwSignal::new(String::new());
    let abund_text = RwSignal::new(String::new());
    let seed_text = RwSignal::new(String::new());

    let read_into = move |sig: RwSignal<String>| {
        move |ev: web_sys::Event| {
            let input: HtmlInputElement = event_target(&ev);
            let Some(file) = input.files().and_then(|f| f.get(0)) else {
                return;
            };
            let gfile = gloo_file::File::from(file);
            spawn_local(async move {
                if let Ok(text) = gloo_file::futures::read_as_text(&gfile).await {
                    sig.set(text);
                }
            });
        }
    };

    let load_custom = move |_| {
        let (c, a, s) = (contigs_text.get(), abund_text.get(), seed_text.get());
        if c.is_empty() || a.is_empty() {
            status.set("Need at least contigs + abundance.".to_string());
            return;
        }
        let ds = Dataset::from_inputs(&c, &a, &s);
        status.set(format!("custom: {} contigs", ds.contigs.len()));
        dataset.set(Some(ds));
    };

    view! {
        <main style="font:15px/1.5 system-ui,sans-serif; max-width:44rem; margin:2.5rem auto; padding:0 1rem;">
            <h1 style="font-size:1.3rem; margin-bottom:0.1rem;">"maxbin-rs · EM inputs"</h1>
            <p style="margin-top:0; opacity:0.7;">
                "Contigs plotted by abundance and length. Seed contigs (initial bins) in orange."
            </p>

            <div style="display:flex; gap:0.75rem; align-items:center; margin:1rem 0;">
                <label>
                    "dataset: "
                    <select on:change=move |ev| {
                        let v = event_target_value(&ev);
                        if v == "custom" {
                            custom.set(true);
                            status.set("Upload contigs + abundance (+ optional seed).".to_string());
                        } else {
                            custom.set(false);
                            load_builtin("bfragilis");
                        }
                    }>
                        <option value="">"— choose —"</option>
                        <option value="bfragilis">"B. fragilis (single genome, smoke test)"</option>
                        <option value="custom">"Custom (upload)"</option>
                    </select>
                </label>
                <span style="opacity:0.7;">{move || status.get()}</span>
            </div>

            <Show when=move || custom.get()>
                <fieldset style="margin-bottom:1rem; border:1px solid currentColor; border-radius:6px; padding:0.75rem;">
                    <legend style="opacity:0.7;">"EM inputs (read in-page, never uploaded)"</legend>
                    <div style="display:grid; grid-template-columns:auto 1fr; gap:0.4rem 0.6rem; align-items:center;">
                        <span>"contigs (FASTA)"</span>
                        <input type="file" on:change=read_into(contigs_text) />
                        <span>"abundance"</span>
                        <input type="file" on:change=read_into(abund_text) />
                        <span>"seed (optional)"</span>
                        <input type="file" on:change=read_into(seed_text) />
                    </div>
                    <button on:click=load_custom style="margin-top:0.6rem;">"Load"</button>
                </fieldset>
            </Show>

            <canvas
                node_ref=canvas_ref
                width=CANVAS_W
                height=CANVAS_H
                style="display:block; width:100%; border:1px solid currentColor; border-radius:6px;"
            />
        </main>
    }
}

/// Fetch the three EM-input files for a built-in dataset under `base/`.
async fn fetch_triple(base: &str) -> Result<(String, String, String), String> {
    async fn get(url: String) -> Result<String, String> {
        gloo_net::http::Request::get(&url)
            .send()
            .await
            .map_err(|e| e.to_string())?
            .text()
            .await
            .map_err(|e| e.to_string())
    }
    let contigs = get(format!("{base}/contigs.fa")).await?;
    let abund = get(format!("{base}/abund")).await?;
    // Seed is optional; treat a fetch failure as "no seeds".
    let seed = get(format!("{base}/seed")).await.unwrap_or_default();
    Ok((contigs, abund, seed))
}

/// Scatter plot: x = abundance (log), y = length (log). Seeds in orange.
fn draw(ctx: &CanvasRenderingContext2d, ds: Option<&Dataset>) {
    let (w, h) = (CANVAS_W as f64, CANVAS_H as f64);
    ctx.set_fill_style_str("#0d1117");
    ctx.fill_rect(0.0, 0.0, w, h);

    let Some(ds) = ds else {
        ctx.set_fill_style_str("#8b949e");
        ctx.set_font("14px system-ui, sans-serif");
        ctx.set_text_align("center");
        let _ = ctx.fill_text("no dataset loaded", w / 2.0, h / 2.0);
        return;
    };

    // Plot only contigs with positive abundance and length (log axes).
    let pts: Vec<&Contig> = ds
        .contigs
        .iter()
        .filter(|c| c.abundance > 0.0 && c.length > 0)
        .collect();
    if pts.is_empty() {
        return;
    }

    let pad = 48.0;
    let log = |v: f64| v.max(1e-9).ln();
    let (mut ax0, mut ax1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ly0, mut ly1) = (f64::INFINITY, f64::NEG_INFINITY);
    for c in &pts {
        ax0 = ax0.min(log(c.abundance));
        ax1 = ax1.max(log(c.abundance));
        ly0 = ly0.min(log(c.length as f64));
        ly1 = ly1.max(log(c.length as f64));
    }
    let sx = |v: f64| pad + (log(v) - ax0) / (ax1 - ax0).max(1e-9) * (w - 2.0 * pad);
    let sy = |v: f64| (h - pad) - (log(v) - ly0) / (ly1 - ly0).max(1e-9) * (h - 2.0 * pad);

    // Axes.
    ctx.set_stroke_style_str("#30363d");
    ctx.set_line_width(1.0);
    ctx.begin_path();
    ctx.move_to(pad, pad);
    ctx.line_to(pad, h - pad);
    ctx.line_to(w - pad, h - pad);
    ctx.stroke();
    ctx.set_fill_style_str("#8b949e");
    ctx.set_font("12px system-ui, sans-serif");
    ctx.set_text_align("center");
    let _ = ctx.fill_text("abundance (log) →", w / 2.0, h - 12.0);
    let _ = ctx.fill_text("length (log) ↑", pad, 18.0);

    // Points: seeds last so they sit on top.
    let mut draw_pt = |c: &Contig, fill: &str, r: f64| {
        ctx.set_fill_style_str(fill);
        ctx.begin_path();
        let _ = ctx.arc(sx(c.abundance), sy(c.length as f64), r, 0.0, TAU);
        ctx.fill();
    };
    for c in pts.iter().filter(|c| !c.is_seed) {
        draw_pt(c, "rgba(79,140,255,0.55)", 3.0);
    }
    for c in pts.iter().filter(|c| c.is_seed) {
        draw_pt(c, "#ff9f40", 5.0);
    }
}
