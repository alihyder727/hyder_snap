restore, 'cirs2000.sav'
openw, lun, 'cirs2000.txt', /get_lun
printf, lun, '# Pressure', n_elements(press)
printf, lun, press[*]
printf, lun, '# Latitude', n_elements(cirs_lat)
printf, lun, cirs_lat[*]
printf, lun, '# Temperature profile', n_elements(cirs_temp)
printf, lun, cirs_temp[*, 0, *]
printf, lun, '# NH3 profile', n_elements(cirs_nh3profile)
printf, lun, cirs_nh3profile[*, 1, *]
Free_lun, lun

end
