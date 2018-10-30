# Copyright 1999-2018 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2

EAPI=6

inherit cmake-utils

DESCRIPTION="Very simple FITS visualisation tool"
HOMEPAGE="https://github.com/eddyem/fits_filter"
SRC_URI=""

LICENSE="GPL-3"
SLOT="0"
KEYWORDS="~amd64 ~x86"
IUSE="cuda leptonica gsl"

DEPEND="
    virtual/opengl
    x11-libs/gtkglext
    sci-libs/cfitsio
    sci-libs/fftw
    cuda? ( dev-util/nvidia-cuda-toolkit )
    leptonica? ( media-libs/leptonica )
    gsl? ( sci-libs/gsl )
"
RDEPEND="${DEPEND}"
S="${WORKDIR}"

src_unpack() {
    cp -r /home/eddy/C-files/fitsview/* ${S} || die
}

CMAKE_USE_DIR="${S}"
BUILD_DIR=${S} 

src_configure() {
    local mycmakeargs=(
	"-DCMAKE_INSTALL_PREFIX=/usr/local"
	"-DUSE_CUDA=$(usex cuda)"
	"-DUSE_LEPTONICA=$(usex leptonica)"
	"-DUSE_GSL=$(usex gsl)"
    )
     cmake-utils_src_configure
}
