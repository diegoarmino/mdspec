./mdspec input tseries.dat spectra.dat
ruby -e'puts readlines.map(&:split).transpose.map{|x|x*" "}' spectra.dat > tspectra.dat
