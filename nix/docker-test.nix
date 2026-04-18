# NixOS VM test — verify the Docker image works end-to-end.
#
# Boots a VM with Docker, loads the Nix-built image, runs maxbin-rs
# on the B. fragilis test dataset, and checks that output bins are produced.
#
#   nix build .#dockerTest

{
  pkgs,
  dockerImage,
  datasets,
}:

pkgs.testers.runNixOSTest {
  name = "maxbin-rs-docker";

  nodes.machine =
    { ... }:
    {
      virtualisation = {
        docker.enable = true;
        diskSize = 4096;
        memorySize = 2048;
      };
    };

  testScript = ''
    machine.wait_for_unit("docker.service")

    # Load the Nix-built image
    machine.succeed("docker load < ${dockerImage}")
    machine.succeed("docker images | grep maxbin-rs")

    # Prepare test data in a writable directory
    machine.succeed("mkdir -p /tmp/testdata /tmp/output")
    machine.succeed("cp ${datasets.bfragilis.contigs} /tmp/testdata/contigs.fa.gz")
    machine.succeed("cp ${datasets.bfragilis.reads1} /tmp/testdata/reads1.fastq.gz")

    # Smoke test: --version exits cleanly
    machine.succeed(
        "docker run --rm ghcr.io/werner291/maxbin-rs:latest --version"
    )

    # Integration test: run full pipeline on B. fragilis
    result = machine.succeed(
        "docker run --rm "
        "-v /tmp/testdata:/input:ro "
        "-v /tmp/output:/output "
        "ghcr.io/werner291/maxbin-rs:latest "
        "--contig /input/contigs.fa.gz "
        "--reads /input/reads1.fastq.gz "
        "--out /output/bins "
        "--thread 1 "
        "2>&1 || true"
    )
    print(result)

    # Verify output bins were created
    machine.succeed("ls /tmp/output/bins.*.fasta")
    machine.succeed("test -f /tmp/output/bins.summary")
  '';
}
