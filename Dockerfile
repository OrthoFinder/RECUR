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
COPY --from=builder /usr/src/recur .


# Create a non-root user, create ExampleData directory, and assign ownership
RUN useradd -m -s /bin/bash recur_user
# RUN useradd -m -s /bin/bash recur_user && \
# mkdir -p /usr/src/recur/ExampleData && \
# chown -R recur_user:recur_user /usr/src/recur

USER recur_user

ENV PATH="/usr/src/recur/src/recur/bin:${PATH}"
EXPOSE 80
ENV NAME=recur

ENTRYPOINT ["python3", "recur.py"]
CMD ["--help"]

# docker container run -it --rm -v $(pwd)/ExampleData:/usr/src/RECUR/ExampleData:Z orthofinder/recur:v1.0.0 -f ExampleData/example_alignments.aln -st AA --outgroups ExampleData/example_alignments.outgroups.txt
