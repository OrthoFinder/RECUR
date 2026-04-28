FROM python:3.13-slim-bookworm AS builder

WORKDIR /usr/src/recur

ENV PIP_DEFAULT_TIMEOUT=120 \
    PIP_RETRIES=10 \
    PIP_NO_CACHE_DIR=1

COPY requirements.txt .

RUN pip install \
    --index-url https://pypi.org/simple \
    -r requirements.txt

COPY . .

RUN chmod +x ./src/recur/bin/iqtree3 || true

RUN cp -r ExampleData /usr/src/recur/default_exampledata


FROM python:3.13-slim-bookworm

WORKDIR /usr/src/recur

# Copy installed packages and code
COPY --from=builder /usr/local/lib/python3.13/site-packages /usr/local/lib/python3.13/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/src/recur /usr/src/recur

# Install gosu
RUN apt-get update \
 && apt-get install -y --no-install-recommends gosu \
 && rm -rf /var/lib/apt/lists/*

RUN chmod -R a+rX /usr/src/recur

RUN echo '#!/bin/bash\n\
set -e\n\
\n\
if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    if ! id -u hostuser &>/dev/null; then\n\
        useradd -o -u "$LOCAL_UID" -m hostuser 2>/dev/null || true\n\
        groupmod -g "$LOCAL_GID" hostuser 2>/dev/null || true\n\
    fi\n\
fi\n\
\n\
if [ -z "$(ls -A ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData} 2>/dev/null)" ]; then\n\
    echo "Populating ExampleData from default_exampledata..." >&2\n\
    cp -r /usr/src/recur/default_exampledata/* ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData}\n\
fi\n\
\n\
if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    echo "chowning /usr/src/recur to $LOCAL_UID:$LOCAL_GID..." >&2\n\
    chown -R $LOCAL_UID:$LOCAL_GID /usr/src/recur || true\n\
fi\n\
\n\
echo "Running as user: $(id)" >&2\n\
\n\
if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    exec gosu $LOCAL_UID:$LOCAL_GID python3 /usr/src/recur/recur.py "$@"\n\
else\n\
    exec python3 /usr/src/recur/recur.py "$@"\n\
fi\n\
' > /usr/local/bin/entrypoint.sh

RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"