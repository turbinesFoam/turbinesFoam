#!/bin/bash
# -----------------------------------------------------------------------------
# Add openfoam debian/ubuntu repository
#
# Copyright (C) 2020-2023 OpenCFD Ltd.
# SPDX-License-Identifier: (GPL-3.0+)
#
# Usage
#     curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
#     wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | sudo bash
# -----------------------------------------------------------------------------

# Input locations
pubkey_url="https://dl.openfoam.com/pubkey.gpg"
repos_url="https://dl.openfoam.com/repos/deb"

# Output names for trusted.gpg.d and sources.list.d
apt_pubkey_name="openfoam.gpg"
apt_source_name="openfoam.list"


# Emit cannot install XXX information and exit
fatal_cannot_install()
{
    echo "Unable to install: $@"
    echo "Your base system has a problem - this should not happen."
    echo "Check your default package repositories."
    echo "Repository installation aborted."
    exit 1
}

# Retrieve the architecture name (amd64, arm64, ...)
# ---------------------------------------------------
unset archName

archName="$(dpkg --print-architecture 2>/dev/null)"
: "${archName:=amd64}"  # Hard-coded failsafe value


# Require code-name
# -----------------
unset codeName

if [ -e /etc/os-release ]
then
    # Check for Ubuntu (or derivatives) first
    [ -n "$codeName" ] || codeName="$(sed -ne 's/^UBUNTU_CODENAME=//p' /etc/os-release)"
    [ -n "$codeName" ] || codeName="$(sed -ne 's/^VERSION_CODENAME=//p' /etc/os-release)"
fi
if [ -z "$codeName" ]
then
    if [ -e /etc/lsb-release ]
    then
        codeName="$(sed -ne 's/^DISTRIB_CODENAME=//p' /etc/lsb-release)"
    fi
    # Not sure why this might work if all others failed, but check anyhow
    [ -n "$codeName" ] || codeName="$(lsb_release -cs 2>/dev/null)"
fi

case "$codeName" in
('' | "n/a")
    codeName=stable
    echo "Unknown distribution code-name, assuming '$codeName'"
    ;;
(*)
    echo "Detected distribution code-name: $codeName"
    ;;
esac

# -----------------


# Need gpg
if ! command -v gpg >/dev/null
then
    echo "Installing gnupg for handling repository keys..."
    apt-get install -y gnupg || fatal_cannot_install gpg
fi


# Need curl or wget to fetch the pubkey.
# If neither are available, install wget (smaller dependency footprint)

unset have_curl have_wget

if command -v curl >/dev/null
then
    have_curl=true
elif command -v wget >/dev/null
then
    have_wget=true
else
    echo "Installing wget..."
    apt-get install -q -y wget
    if command -v wget >/dev/null
    then
        have_wget=true
    else
        fatal_cannot_install wget
    fi
fi


# pubkey registry
# ---------------
# eval $(apt-config shell APT_TRUSTED_PARTS Dir::Etc::trustedparts/d)
# apt_pubkey_path="${APT_TRUSTED_PARTS}${apt_pubkey_name}"
apt_pubkey_path="/etc/apt/trusted.gpg.d/${apt_pubkey_name}"

# sources list registry
# ---------------------
# eval $(apt-config shell APT_SOURCE_PARTS Dir::Etc::sourceparts/d)
# apt_source_path="${APT_SOURCE_PARTS}${apt_source_name}"
apt_source_path="/etc/apt/sources.list.d/${apt_source_name}"


fileAction="Added"
if [ -f "$apt_source_path" ]
then
    fileAction="Overwrote"
fi

echo "### THIS FILE IS AUTOMATICALLY CONFIGURED ###
# You may comment out this entry, but any other modifications may be lost.
deb [arch=${archName}] ${repos_url} ${codeName} main" \
    > "$apt_source_path"
echo "$fileAction $apt_source_path"


fileAction="Added"
if [ -f "$apt_pubkey_path" ]
then
    fileAction="Overwrote"
fi

echo -n "Importing openfoam gpg key... "
if [ -n "$have_curl" ]
then
    curl -L "${pubkey_url}" 2>/dev/null | gpg --dearmor > "$apt_pubkey_path"
else
    wget -O - -nv "${pubkey_url}" 2>/dev/null | gpg --dearmor > "$apt_pubkey_path"
fi

if [ -s "$apt_pubkey_path" ]
then
    echo "done"
    echo "$fileAction $apt_pubkey_path"
else
    echo "Error"
    echo "Failed to add $apt_pubkey_path"
    exit 1
fi


# Update
echo -n "Running apt-get update... "
apt-get update >/dev/null 2>/dev/null
echo "done"

echo
echo "The repository is setup! You can now install packages."

# ---------------------------------------------------------------------------
