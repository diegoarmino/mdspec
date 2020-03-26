ruby -e'puts readlines.map(&:split).transpose.map{|x|x*" "}' spectra.dat > transposed-spectra.dat
