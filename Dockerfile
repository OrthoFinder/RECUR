FROM python:3.12-slim AS builder

WORKDIR /usr/src/recur

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN chmod +x ./src/recur/bin/iqtree2 || true

RUN cp -r ExampleData /usr/src/recur/default_exampledata

FROM python:3.12-slim

WORKDIR /usr/src/recur

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/src/recur /usr/src/recur

RUN apt-get update && apt-get install -y gosu && rm -rf /var/lib/apt/lists/*

RUN chmod -R a+rX /usr/src/recur

#
# ENTRYPOINT that:
#  - runs as root
#  - copies ExampleData if empty
#  - chowns everything to $LOCAL_UID:$LOCAL_GID if provided
#  - finally 'gosu' to that user & run python3
#
RUN echo '#!/bin/bash\n\
    set -e\n\
    \n\
    # 1) If user specified LOCAL_UID and LOCAL_GID, create that user if needed.\n\
    if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    # Create user if it does not exist\n\
    if ! id -u hostuser &>/dev/null; then\n\
    useradd -o -u "$LOCAL_UID" -m hostuser 2>/dev/null || true\n\
    groupmod -g "$LOCAL_GID" hostuser 2>/dev/null || true\n\
    fi\n\
    fi\n\
    \n\
    # 2) If the ExampleData folder is empty, copy from default_exampledata.\n\
    if [ -z "$(ls -A ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData} 2>/dev/null)" ]; then\n\
    echo "Populating ExampleData from default_exampledata..." >&2\n\
    cp -r /usr/src/recur/default_exampledata/* ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData}\n\
    fi\n\
    \n\
    # 3) Chown the code + data so our host user can remove/modify it on the host.\n\
    if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    echo "chowning /usr/src/recur to $LOCAL_UID:$LOCAL_GID..." >&2\n\
    chown -R $LOCAL_UID:$LOCAL_GID /usr/src/recur || true\n\
    fi\n\
    \n\
    echo "Running as user: $(id)" >&2\n\
    \n\
    # 4) Finally drop privileges if LOCAL_UID:GID set, else just run as root.\n\
    if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    exec gosu $LOCAL_UID:$LOCAL_GID python3 /usr/src/recur/recur.py "$@"\n\
    else\n\
    exec python3 /usr/src/recur/recur.py "$@"\n\
    fi\n\
    ' > /usr/local/bin/entrypoint.sh

RUN chmod +x /usr/local/bin/entrypoint.sh

# IMPORTANT: do NOT override with 'USER 1000:1000' or similar
# We want to start as root so we can create users/chown. Then drop privileges with gosu.
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]




# docker rmi -f $(docker images -aq)
# UID=$(id -u) GID=$(id -g) docker compose up --build
#  docker tag orthofinder/recur:v1.0.0 orthofinder/recur:latest
#  docker push orthofinder/recur:v1.0.0
#  docker push orthofinder/recur:latest
# docker run -it --rm --entrypoint bash orthofinder/recur:v1.0.0

# docker container run -it --rm orthofinder/recur:v1.0.0
# docker run -it --rm \
#     -v $(pwd)/MyData:/usr/src/recur/MyData \
#     -e LOCAL_UID=$(id -u) -e LOCAL_GID=$(id -g) \
#     orthofinder/recur:v1.0.0 \
#     -f MyData/example_alignments.aln \
#     -st AA \
#     --outgroups MyData/example_alignments.outgroups.txt

