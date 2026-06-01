//! Human-facing prose, in one place.
//!
//! Everything in this file is meant to be written in Werner's voice, not
//! Claude's. The strings below are DRAFT scaffolding only - placeholders that
//! say roughly the right thing so the layout has text to render. Rewrite them.
//! Nothing here is reviewed copy until it has been rewritten by hand.
//!
//! The rest of the crate (parsing in `data.rs`, layout/metadata in `main.rs`)
//! is machinery; the words live here.

/// Eyebrow line above the title.
pub const KICKER: &str = "maxbin-rs · the binning algorithm, on real data";

/// Masthead title (the project name).
pub const PAGE_TITLE: &str = "maxbin-rs";

/// Opening, written by Werner. Spelling fixed (MaxBin2, succeeded), split into
/// two paragraphs for readability; `*understand*` renders as emphasis. Words
/// are otherwise his.
pub const INTRO_1: &str = "A few weeks ago, I implemented maxbin-rs, a full \
    Rust rewrite of MaxBin2. MaxBin2 is a metagenomics binning tool, which \
    takes a set of contigs, and classifies them into coherent \"bins\", i.e. \
    genomes belonging to the same or similar species. At the time, I treated \
    it entirely abstractly - I got byte-identical inputs and outputs, trimmed \
    dependencies, went for an exact rewrite. It was meant to be the same tool, \
    which succeeded. It was faster (mostly due to using 64 bit floats and \
    taking advantage of LLVM's automatic inlining), had fewer flaky \
    dependencies that were properly locked in place.";

pub const INTRO_2: &str = "However, what bothered me to this day is that I \
    didn't quite fully *understand* this tool, and I cannot seriously pretend \
    to be its maintainer if I do not know it inside and out. Therefore, I \
    decided to visualize its internals, the result of which can be seen here.";

/// Section heading.
pub const TITLE: &str = "1 · What goes in";

/// Section-1 opening, written by Werner. Spelling fixed (pipeline, rewrote,
/// although, reimplemented, FragGeneScanRs); the parenthetical "(pipeline
/// stages)" placeholders kept his words but dropped the stray parens; split at
/// his ` - ` breaks into paragraphs. Words otherwise his.
pub const LEDE_1: &str = "The original MaxBin2 was a pipeline that took \
    assembled DNA fragments (contigs) and sorted them into \"bins\" \
    (groupings) based on similarity, ideally one bin per species. It consisted \
    of two things: a set of pipeline stages, and a custom \
    expectation-maximization algorithm with a custom distance metric.";

pub const LEDE_2: &str = "The EM part is most of what we rewrote here, as well \
    as the Perl-based orchestration layer that actually ran the pipeline \
    including the other stages. Those original stages were not rewritten, \
    although I did replace FragGeneScan with FragGeneScanRs, a Rust rewrite of \
    that tool.";

pub const LEDE_3: &str = "The page you see here lets you run the parts we \
    reimplemented in your browser. It's the real algorithm; no part of this is \
    pre-rendered, so you can run it with any dataset. First, pick one below.";

/// Heading for the pipeline-overview diagram.
pub const PIPELINE_HEADING: &str = "The pipeline, and where this dataset enters";

/// Intro above the pipeline diagram. DRAFT, rewrite in your hand.
pub const PIPELINE_INTRO: &str = "DRAFT, rewrite in your hand: \
    maxbin is a chain of steps. Only the last one, the EM, runs in your \
    browser; the earlier steps are external tools (a read mapper, a gene \
    caller, HMMER) that produce files this page reads but cannot run. Each \
    dataset ships some of those files already made, so it enters the chain \
    partway down. Below: the full chain, what each step does, and which of \
    its outputs this dataset hands you.";

/// Blurb shown under the B. fragilis dataset.
pub const BFRAGILIS_BLURB: &str = "DRAFT, rewrite in your hand: \
    The nf-core/modules MaxBin2 test sample: a single bacterial genome, \
    sequenced once. Small enough that the whole algorithm can later run live \
    in this page.";

/// Blurb shown under the CAMI I High dataset.
pub const CAMI_BLURB: &str = "DRAFT, rewrite in your hand: \
    The standard metagenome binning benchmark (Sczyrba et al., Nature Methods \
    2017): a synthetic community sequenced across five samples. Counts and \
    lengths below are read straight from the depth file.";

/// Note shown for CAMI, explaining why the assembly FASTA is not bundled.
pub const CAMI_NOTE: &str = "DRAFT, rewrite in your hand: \
    The 2.80 Gbp assembly FASTA (~817 MB gzipped) is the bulk input and is not \
    bundled here; the depth file alone is enough to describe the contig set. \
    The full algorithm will run in-browser on the smaller dataset.";

/// Explanatory clause appended after the seed count.
pub const SEED_DESC: &str = "the initial group centers the algorithm starts from";

/// Shown when a dataset has no bundled seed file.
pub const NO_SEED_NOTE: &str = "DRAFT, rewrite in your hand: \
    No seed file bundled for this dataset; seeds come from the marker-gene \
    step (HMMER + gene caller), which does not run in the browser.";
