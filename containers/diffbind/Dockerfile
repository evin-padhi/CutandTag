FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/diffbind.yml /tmp/diffbind.yml
RUN micromamba install --yes --name base --file /tmp/diffbind.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER bin/diffbind_analysis.R /opt/nanocut/bin/diffbind_analysis.R
ENV PATH="/opt/conda/bin:${PATH}"
ENTRYPOINT []
CMD ["Rscript", "/opt/nanocut/bin/diffbind_analysis.R"]
