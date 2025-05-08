FROM mambaorg/micromamba:2.1.1
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER kmer_filter.py /opt/kmer_filter.py
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/kmer_filter.py"]
