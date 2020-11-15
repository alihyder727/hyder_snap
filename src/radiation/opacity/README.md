# O2-O2_2011.cia, O2-O2_2011_cheng.cia
The original file is download from
https://www.cfa.harvard.edu/HITRAN/HITRAN2012/CIA/Main-Folder/O2-O2/

It provides 3 data sets between 7450 - 8487 wavenumbers. However, those three
data sets have three different wavenumber grids. The first one starts from
7450.38 and ends at 8477.18 with 4194 wavenumbers in total. The second
one starts from 7500.089 and ends at 8486.485 with 4029 wavenumbers in total.
The third one starts from 7450.132 and ends at 8487.465 with 4237 wavenumbers in
total. Having three slightly different wavenumber grids is hard to manipulate in
the program. To fix this, Cheng Li deleted the last two data fields in this band
and renamed the data file to O2-O2_2011_cheng.cia .
