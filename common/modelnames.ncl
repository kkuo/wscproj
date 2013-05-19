modelfnames = (/"CCSM4",\
				"CNRM-CM5",\
				"GFDL-ESM2G", \
				"INMCM4",\
				"IPSL-CM5A-LR",\
				"HadESM2ES",\
				"MIROC5",\
				"MIROC-ESM",\
				"MIROC-ESM-CHEM",\
				"MPI-ESM-LR",\
				"MRI-CGCM3",\
				"NorESM"/)

husfnames = (/	"CCSM4/keep/hus_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
				"CNRM-CM5/hus_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
				"GFDL-ESM2G/hus.nc", \
				"INMCM4/hus_Amon_inmcm4_1970_2005.nc",\
				"IPSL-CM5A-LR/hus_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
				"HadESM2ES/hus_195912-200511.nc",\
				"MIROC5/hus_Amon_MIROC5_1970_2012.nc",\
				"MIROC-ESM/hus_Amon_MIROC-ESM_1850_2005.nc",\
				"MIROC-ESM-CHEM/hus_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
				"MPI-ESM-LR/hus_1970_to_2005.nc",\
				"MRI-CGCM3/hus_Amon_197001-200512.nc",\
				"NorESM/hus_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

uafnames = (/	"CCSM4/keep/ua_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
				"CNRM-CM5/ua_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
				"GFDL-ESM2G/ua.nc", \
				"INMCM4/ua_Amon_inmcm4_1970_2005.nc",\
				"IPSL-CM5A-LR/ua_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
				"HadESM2ES/ua_195912-200511.nc",\
				"MIROC5/ua_Amon_MIROC5_1970_2012.nc",\
				"MIROC-ESM/ua_Amon_MIROC-ESM_1850_2005.nc",\
				"MIROC-ESM-CHEM/ua_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
				"MPI-ESM-LR/ua_1970_to_2005.nc",\
				"MRI-CGCM3/ua_Amon_197001-200512.nc",\
				"NorESM/ua_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

vafnames = (/	"CCSM4/keep/va_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
				"CNRM-CM5/va_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
				"GFDL-ESM2G/va.nc", \
				"INMCM4/va_Amon_inmcm4_1970_2005.nc",\
				"IPSL-CM5A-LR/va_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
				"HadESM2ES/va_195912-200511.nc",\
				"MIROC5/va_Amon_MIROC5_1970_2012.nc",\
				"MIROC-ESM/va_Amon_MIROC-ESM_1850_2005.nc",\
				"MIROC-ESM-CHEM/va_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
				"MPI-ESM-LR/va_1970_to_2005.nc",\
				"MRI-CGCM3/va_Amon_197001-200512.nc",\
				"NorESM/va_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

undef("getmodelcount")
function getmodelcount()
begin
	return dimsizes(husfnames)
end

undef("gethusmodelbyindex")
function gethusmodelbyindex(idx:numeric)
begin
	size = getmodelcount()
	if(idx .ge. size) then
		return ""
	end if
	return husfnames(idx)
end

undef("getuamodelbyindex")
function getuamodelbyindex(idx:numeric)
begin
	size = getmodelcount()
	if(idx .ge. size) then
		return ""
	end if
	return uafnames(idx)
end

undef("getvamodelbyindex")
function getvamodelbyindex(idx:numeric)
begin
	size = getmodelcount()
	if(idx .ge. size) then
		return ""
	end if
	return vafnames(idx)
end

undef("getmodelnamebyindex")
function getmodelnamebyindex(idx:numeric)
begin
	size = getmodelcount()
	if(idx .ge. size) then
		return ""
	end if
	return modelfnames(idx)
end

