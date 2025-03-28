# Builder stage
FROM python:3.12-slim AS builder

WORKDIR /usr/src/recur

# Copy only the requirements first to leverage Docker cache
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Now copy the rest of the application files
COPY . .

# Ensure iqtree2 is executable if it exists
RUN chmod +x ./src/recur/bin/iqtree2 || true

# Final stage
FROM python:3.12-slim

WORKDIR /usr/src/recur

# Copy installed dependencies and binaries from the builder stage
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy application source code
COPY --from=builder /usr/src/recur .

# Create a non-root user, create ExampleData directory, and assign ownership
RUN useradd -m -s /bin/bash recur_user && \
    mkdir -p /usr/src/recur/ExampleData && \
    chown -R recur_user:recur_user /usr/src/recur

USER recur_user

# Update PATH to include the bin directory
ENV PATH="/usr/src/recur/src/recur/bin:${PATH}"

# Expose port 80
EXPOSE 80

# Define environment variable
ENV NAME=recur

# Run recur.py when the container launches
CMD ["python3", "./recur.py"]
