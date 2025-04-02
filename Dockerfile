# -------------------------------------------
# 1) Builder stage: install dependencies
# -------------------------------------------
FROM python:3.12-slim AS builder

WORKDIR /usr/src/recur

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy all source code (including ExampleData/) to build context
COPY . .

# For demonstration, if you have any binaries you want to chmod:
RUN chmod +x ./src/recur/bin/iqtree2 || true

# Copy the "shipped" ExampleData into a safe location
RUN cp -r ExampleData /usr/src/recur/default_exampledata


# -------------------------------------------
# 2) Final stage
# -------------------------------------------
FROM python:3.12-slim
WORKDIR /usr/src/recur

# Copy Python site-packages and local bin from builder
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy the entire /usr/src/recur from builder
COPY --from=builder /usr/src/recur /usr/src/recur

# Install gosu for UID/GID switching
RUN apt-get update && apt-get install -y gosu && rm -rf /var/lib/apt/lists/*

# Make everything readable/executable as needed
RUN chmod -R a+rX /usr/src/recur

# -------------------------------------------
# Entry point script:
# - if RECUR_DATA_DIR is empty, copy default_exampledata there
# - optionally switch to LOCAL_UID:LOCAL_GID
# - run python3 recur.py with passed args
# -------------------------------------------
RUN echo '#!/bin/bash\n\
    set -e\n\
    \n\
    # Optional: if user supplies LOCAL_UID and LOCAL_GID, switch to that user.\n\
    if [ -n "$LOCAL_UID" ] && [ -n "$LOCAL_GID" ]; then\n\
    # Create a user "user" if it does not exist\n\
    if ! id -u user &>/dev/null; then\n\
    useradd -u "$LOCAL_UID" -o -m user 2>/dev/null || true\n\
    groupmod -g "$LOCAL_GID" user 2>/dev/null || true\n\
    fi\n\
    # Make sure /usr/src/recur is owned by this user\n\
    chown -R user:user /usr/src/recur || true\n\
    echo "Switching to UID:GID $LOCAL_UID:$LOCAL_GID..." >&2\n\
    exec gosu user "$0" "$@"\n\
    fi\n\
    \n\
    # If ExampleData dir is empty, copy from default_exampledata.\n\
    if [ -z "$(ls -A ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData} 2>/dev/null)" ]; then\n\
    echo "Populating ExampleData from default_exampledata..." >&2\n\
    cp -r /usr/src/recur/default_exampledata/* ${RECUR_DATA_DIR:-/usr/src/recur/ExampleData}\n\
    fi\n\
    \n\
    echo "Running as user: $(id)" >&2\n\
    exec python3 recur.py "$@"\n' \
    > /usr/local/bin/entrypoint.sh

RUN chmod +x /usr/local/bin/entrypoint.sh

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

