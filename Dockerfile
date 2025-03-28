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

# Set up default non-root user
RUN groupadd -g 1000 recur_group && \
    useradd -u 1000 -g recur_group recur_user && \
    chown -R recur_user:recur_group /usr/src/recur

USER recur_user

ENV PATH="/usr/src/recur/src/recur/bin:${PATH}"
EXPOSE 80
ENV NAME=recur

ENTRYPOINT ["python3", "./recur.py"]
CMD ["--help"]
