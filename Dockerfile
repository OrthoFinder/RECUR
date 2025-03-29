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

RUN apt-get update && apt-get install -y gosu && rm -rf /var/lib/apt/lists/*

RUN chmod -R a+rX /usr/src/recur

RUN echo '#!/bin/bash\n\
    set -e\n\
    \n\
    USER_ID=${LOCAL_UID:-1000}\n\
    GROUP_ID=${LOCAL_GID:-1000}\n\
    \n\
    if ! getent group "$GROUP_ID" > /dev/null 2>&1; then\n\
    groupadd -g "$GROUP_ID" recur_group\n\
    fi\n\
    \n\
    if ! id -u "$USER_ID" > /dev/null 2>&1; then\n\
    useradd -u "$USER_ID" -g "$GROUP_ID" -m recur_user\n\
    fi\n\
    \n\
    RECUR_DATA_DIR=${RECUR_DATA_DIR:-/usr/src/recur/MyData}\n\
    if [ -d "$RECUR_DATA_DIR" ]; then\n\
    chmod -R a+rwX "$RECUR_DATA_DIR" || true\n\
    else\n\
    echo "[RECUR ENTRYPOINT] Skipping chmod, data directory not found: $RECUR_DATA_DIR"\n\
    fi\n\
    \n\
    exec gosu "$USER_ID:$GROUP_ID" python3 recur.py "$@"' \
    > /usr/local/bin/entrypoint.sh && chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]


# docker rmi -f $(docker images -aq)
# UID=$(id -u) GID=$(id -g) docker compose up --build
#  docker tag orthofinder/recur:v1.0.0 orthofinder/recur:latest
#  docker push orthofinder/recur:v1.0.0
#  docker push orthofinder/recur:latest

# docker container run -it --rm orthofinder/recur:v1.0.0
# docker run -it --rm \
#     -v $(pwd)/MyData:/usr/src/recur/MyData \
#     -e LOCAL_UID=$(id -u) -e LOCAL_GID=$(id -g) \
#     orthofinder/recur:v1.0.0 \
#     -f MyData/example_alignments.aln \
#     -st AA \
#     --outgroups MyData/example_alignments.outgroups.txt

