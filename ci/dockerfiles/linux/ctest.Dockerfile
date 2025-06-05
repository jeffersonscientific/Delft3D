FROM containers.deltares.nl/base_linux_containers/8-base:latest

# Install system dependencies and clean up packages afterwards
RUN dnf update --assumeyes && \
    dnf install --assumeyes python3.12 && \
    dnf clean all && \
    rm --recursive --force /var/cache/dnf

RUN dnf install cmake -y
RUN echo "cmake version: $(cmake --version)"
