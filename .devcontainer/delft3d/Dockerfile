FROM containers.deltares.nl/delft3d-dev/delft3d-third-party-libs:oneapi-2024-ifx-release AS base

ARG USER_UID=1000

# Set bash as the default shell for RUN commands.
SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

# Install some tools from the almalinux package manager
RUN --mount=type=cache,target=/var/cache/dnf,sharing=locked <<"EOF"
dnf install -y sudo git dos2unix which less tree strace \
    jq java-17-openjdk java-17-openjdk-devel maven
update-alternatives --set java java-17-openjdk.x86_64
EOF

ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk

# Install Python 3.11.9 from source
RUN <<"EOF"
cd /tmp
curl -O https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz
tar -xzf Python-3.11.9.tgz
cd Python-3.11.9
./configure --enable-optimizations --with-ensurepip=install
make -j$(nproc)
sudo make altinstall
# Clean up
cd /
rm -rf /tmp/Python-3.11.9*
# Create symlinks for system-wide access
sudo ln -sf /usr/local/bin/python3.11 /usr/local/bin/python3
sudo ln -sf /usr/local/bin/python3.11 /usr/local/bin/python
sudo ln -sf /usr/local/bin/pip3.11 /usr/local/bin/pip3
sudo ln -sf /usr/local/bin/pip3.11 /usr/local/bin/pip
EOF

# Create non-root user and grant it sudo privileges:
# See: https://code.visualstudio.com/remote/advancedcontainers/add-nonroot-user
RUN <<"EOF"
useradd --uid $USER_UID -ms /bin/bash dev
echo "dev ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/dev
chmod 0440 /etc/sudoers.d/dev
EOF

# Change user from 'root' to 'dev', and set 'dev's home directory as the image workdir.
USER dev
WORKDIR /home/dev

# Create mount points for bash history and MinIO credentials in the user's home directory.
RUN <<"EOF"
mkdir ~/.history ~/.aws
touch ~/.history/.bash_history
EOF


# Install Python tools using pip
ENV PATH=/usr/local/bin:/home/dev/.local/bin:$PATH
RUN <<"EOF"
python3.11 -m pip install --user fortls fprettify "dvc[s3]"
EOF

# Install git prompt and completion scripts.
ADD --chmod=755 \
    https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash \
    https://raw.githubusercontent.com/git/git/master/contrib/completion/git-prompt.sh \
    /home/dev/.local/bin/

# Configure interactive bash environment.
RUN cat >> /home/dev/.bashrc <<EOF
# Source Intel OneApi 
source /opt/intel/oneapi/setvars.sh

# Set Compilers and library path
export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:\$PKG_CONFIG_PATH
export FC=mpiifx
export CXX=mpicxx
export CC=mpiicx

# Bash history
export PROMPT_COMMAND='history -a'
export HISTFILE=/home/dev/.history/.bash_history

# Git completion and prompt support.
source /home/dev/.local/bin/git-completion.bash
source /home/dev/.local/bin/git-prompt.sh

# Aliases
alias ls='ls --color=auto'

# Set bash prompt.
export PS1='[\u@\h \w\$(__git_ps1 " (%s)")]\n\\$ '
EOF
