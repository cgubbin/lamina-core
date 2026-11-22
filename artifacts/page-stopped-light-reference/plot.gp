$Points << EOD
7.41,805
30.3,778
EOD
set terminal svg size 1000,780 font 'sans,13'
set output 'reference.svg'
set datafile separator ','
set multiplot layout 2,2 title "Page stopped-light waveguide | independent passive TM reference"
set xlabel 'In-plane wavenumber (1/µm)'
set ylabel 'Re ω (rad/ps)'
plot 'dispersion.csv' using 1:2 with lines lw 2 title 'TM branch', $Points using 1:2 with points pt 7 title 'Page, Table 3.2'
set ylabel 'Im ω (rad/ps)'
plot 'dispersion.csv' using 1:3 with lines lw 2 title 'Temporal decay'
set xlabel 'Depth from ITO substrate (µm)'
set ylabel '|H_y|, substrate interface = 1'
set key font ',9'
plot for [b in '6 7.41 9 26 30.3 34'] 'fields.csv' using (column(1)==real(b)?column(2):1/0):3 with lines title sprintf('β = %s',b)
set ylabel '|E_z| / (Z₀ |H_y(0)|)'
plot for [b in '6 7.41 9 26 30.3 34'] 'fields.csv' using (column(1)==real(b)?column(2):1/0):5 with lines title sprintf('β = %s',b)
unset multiplot
