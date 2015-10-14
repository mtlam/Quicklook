# Quicklook

Requires the [PyPulse](https://github.com/mtlam/PyPulse) package to run. For now it is necessary to be in the same directory.

Usage: python quicklook.py [OPTIONS] [PSRFITS filename]

Options
-------
ext - Save Quickllok plot to various image formats (pdf, png, etc.)

min - Plot with time units of minutes [default: False]

nchan - Number of frequency channels to average file down [default: 16]

save - Save the output to a .npz which can later be read in

template - Provide a template file

nodm - do not perform DM calculation

depth - The maximum depth when calculating the DM

iters - Number of DM iterations to perform