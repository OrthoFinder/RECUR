FROM python:3.12-slim AS builder

WORKDIR /usr/src/recur
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
RUN chmod +x ./src/recur/bin/iqtree2 || true

FROM python:3.12-slim
WORKDIR /usr/src/recur

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/src/recur /usr/src/recur

# Create embedded entrypoint script
RUN echo '#!/bin/bash\n\
    set -e\n\
    \n\
    python3 recur.py "$@"\n\
    \n\
    chmod -R a+rwX /usr/src/recur/ExampleData || true' \
    > /usr/local/bin/entrypoint.sh && chmod +x /usr/local/bin/entrypoint.sh

ENV PATH="/usr/src/recur/src/recur/bin:${PATH}"
EXPOSE 80
ENV NAME=recur

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]


# docker container run -it --rm -v $(pwd)/MyData:/usr/src/recur/MyData:Z orthofinder/recur:v1.0.0 -f MyData/example_alignments.aln -st AA --outgroups MyData/example_alignments.outgroups.txt