load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)
	PLEV = 92500

	; define GL region
	; need to be in float for gc_latlon
	slat = 40.
	nlat = 50.
	elon = 285.
	wlon = 265.

	;slat = 41.
	;nlat = 49.
	;elon = 285.;-360
	;wlon = 266.;-360
	PADDING = 20

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12

	;asciidir = "/share/kkuo/AR5/output/"
	;asciifile = "BCC"
	;asciiext = ".txt"

	uafnames = (/	"/share/kkuo/AR5/CCSM4/keep/ua_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/ua_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/ua.nc", \
					"/share/kkuo/AR5/INMCM4/ua_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/ua_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/ua_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/ua_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/ua_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/ua_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/ua_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/ua_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/ua_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

	vafnames = (/	"/share/kkuo/AR5/CCSM4/keep/va_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/va_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/va.nc", \
					"/share/kkuo/AR5/INMCM4/va_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/va_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/va_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/va_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/va_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/va_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/va_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/va_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/va_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

	fu = addfile(uafnames(10),"r")
	fv = addfile(vafnames(10),"r")

	YYYYMM = cd_calendar(fu->time,-1)
	iStart = ind(YYYYMM.eq.(yrStart*100+1))
	iLast = ind(YYYYMM.eq.(yrLast*100+12))
	print(cd_calendar(fu->time(iStart),-1))
	print(cd_calendar(fu->time(iLast),-1))

	uwind1 = fu->ua(time|iStart:iLast,plev|:,lat|:,lon|:)
	uwind = uwind1({time|:},{plev|PLEV},{lat|slat-PADDING:nlat+PADDING},{lon|wlon-PADDING:elon+PADDING})
	delete(uwind1)
	printMinMax(uwind,True)

	vwind1 = fv->va(time|iStart:iLast,plev|:,lat|:,lon|:)
	vwind = vwind1({time|:},{plev|PLEV},{lat|slat-PADDING:nlat+PADDING},{lon|wlon-PADDING:elon+PADDING})
	delete(vwind1)
	printMinMax(vwind,True)

	dims = dimsizes(uwind)
	nmonths = dims(0)
	lat2 = dims(1)
	lon2 = dims(2)
	lats = uwind&lat
	lons = uwind&lon
	print(ismissing(vwind))
	print(ismissing(uwind))
	exit

;************************************************
; calculate the seasonal and annual climatology
;************************************************

	;annual
	nyears = nmonths/12
	uwind_annual = new((/nyears,lat2,lon2/), float)
	vwind_annual = new((/nyears,lat2,lon2/), float)
	month_idx = 0
	do i = 0,nyears-1
		uwind_annual(i,:,:) = dim_avg_n_Wrap(uwind(month_idx:month_idx+11,:,:),0)
		vwind_annual(i,:,:) = dim_avg_n_Wrap(vwind(month_idx:month_idx+11,:,:),0)
		month_idx = month_idx + 12
	end do
	;printVarSummary(uwind_annual)
	;printVarSummary(vwind_annual)

	uwind_lta = dim_avg_n_Wrap(uwind_annual,0)
	vwind_lta = dim_avg_n_Wrap(vwind_annual,0)

	;seasonal
	; winter might be tricky because it is across years

	U_MAM = new((/nyears, lat2,lon2/), float)
	V_MAM = new((/nyears, lat2,lon2/), float)
	U_JJA = new((/nyears, lat2,lon2/), float)
	V_JJA = new((/nyears, lat2,lon2/), float)
	U_SON = new((/nyears, lat2,lon2/), float)
	V_SON = new((/nyears, lat2,lon2/), float)

	iYear = 0
	do i = 2, nmonths-1, 12
		; MAM
		U_MAM(iYear,:,:) = dim_avg_n_Wrap(uwind(i:i+2,:,:),0)
		V_MAM(iYear,:,:) = dim_avg_n_Wrap(vwind(i:i+2,:,:),0)
		; JJA
		U_JJA(iYear,:,:) = dim_avg_n_Wrap(uwind(i+3:i+2+3,:,:),0)
		V_JJA(iYear,:,:) = dim_avg_n_Wrap(vwind(i+3:i+2+3,:,:),0)
		; SON
		U_SON(iYear,:,:) = dim_avg_n_Wrap(uwind(i+6:i+2+6,:,:),0)
		V_SON(iYear,:,:) = dim_avg_n_Wrap(vwind(i+6:i+2+6,:,:),0)
		iYear = iYear + 1
	end do
	; average over all the years
	U_MAM_CLM = dim_avg_n_Wrap(U_MAM,0)
	V_MAM_CLM = dim_avg_n_Wrap(V_MAM,0)
	U_JJA_CLM = dim_avg_n_Wrap(U_JJA,0)
	V_JJA_CLM = dim_avg_n_Wrap(V_JJA,0)
	U_SON_CLM = dim_avg_n_Wrap(U_SON,0)
	V_SON_CLM = dim_avg_n_Wrap(V_SON,0)
	

	;test = new(36, float)
	;test(:) = 0
	;test(0:1) = 1
	;do i = 11, 36-1-12, 12
	;	test(i:i+2) = 1
	;end do
	;test(36-2:35) = 1
	;print(test)

	; the DJF is a little different because the last year only has 1 month,
	; so I am going to not include it in the average
	; the first year is going to be only JF
	U_DJF = new((/nyears+1, lat2,lon2/), float)
	V_DJF = new((/nyears+1, lat2,lon2/), float)

	U_DJF(:,:,:) = 0
	V_DJF(:,:,:) = 0

	iYear = 0
	U_DJF(iYear,:,:) = dim_avg_n_Wrap(uwind(0:1,:,:),0)
	V_DJF(iYear,:,:) = dim_avg_n_Wrap(vwind(0:1,:,:),0)
	iYear = 1
	
	do i = 11, nmonths-1-12, 12
		; DJF
		U_DJF(iYear,:,:) = dim_avg_n_Wrap(uwind(i:i+2,:,:),0)
		V_DJF(iYear,:,:) = dim_avg_n_Wrap(vwind(i:i+2,:,:),0)
		iYear = iYear + 1
	end do
	printVarSummary(U_DJF)

	U_DJF_CLM = dim_avg_n_Wrap(U_DJF,0)
	V_DJF_CLM = dim_avg_n_Wrap(V_DJF,0)
	
	
;************************************************
; plotting parameters
;************************************************

; add polygons to the plot
	xpts = (/wlon, wlon, elon, elon, wlon/)
	ypts = (/slat, nlat, nlat, slat, slat/)

	respoly = True
	respoly@gsLineColor = "red"
	dum = new(40,graphic)


	;sPLEV = sprinti("%i",PLEV)
	;strs = (/"UV_",sPLEV/)
	;fname = str_concat(strs)
	wks   = gsn_open_wks ("pdf","MRI")                ; open workstation
	gsn_define_colormap(wks,"gui_default")


	vcres                         = True             ; plot mods desired

	vcres@gsnFrame                = False            ; so we can draw time stamp
	vcres@gsnDraw                 = False            ; so we can draw time stamp
	vcres@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
	vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
	vcres@vcRefLengthF            = 0.045            ; define length of vec ref
	vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
	vcres@vcMinDistanceF          = 0.017
	vcres@gsnAddCyclic = False
	vcres@gsnMaximize = False
	vcres@gsnSpreadColors = True
	vcres@vcMonoLineArrowColor = True
	vcres@mpMinLonF = wlon-PADDING
	vcres@mpMaxLonF = elon+PADDING
	vcres@mpMinLatF = slat-PADDING
	vcres@mpMaxLatF = nlat+PADDING
	vcres@tiMainString = "Annual (92500)"
	;vcres@cnLinesOn = True
	;vcres@cnFillOn = True
	vcres@mpFillOn = False
	;vcres@mpGeophysicalLineThicknessF = 1.5


	;uwind_lta@long_name = ""
	;vwind_lta@long_name = ""
	plot1 = new(5,graphic)

	vcres@gsnLeftString = ""
	plot1(0) = gsn_csm_vector_map_ce(wks,uwind_lta,vwind_lta,vcres)
	vcres@tiMainString = "Spring"
	plot1(1) = gsn_csm_vector_map_ce(wks,U_MAM_CLM,V_MAM_CLM,vcres)
	vcres@tiMainString = "Summer"
	plot1(2) = gsn_csm_vector_map_ce(wks,U_JJA_CLM,V_JJA_CLM,vcres)
	vcres@tiMainString = "Autum"
	plot1(3) = gsn_csm_vector_map_ce(wks,U_SON_CLM,V_SON_CLM,vcres)
	vcres@tiMainString = "Winter"
	plot1(4) = gsn_csm_vector_map_ce(wks,U_DJF_CLM,V_DJF_CLM,vcres)

	;plot(0) = gsn_csm_vector_map_ce(wks,uwind_lta,vwind_lta,vcres)

	delete(uwind)
	delete(vwind)

	; now 85000 level
	PLEV = 85000
	uwind1 = fu->ua(time|iStart:iLast,plev|:,lat|:,lon|:)
	uwind = uwind1({time|:},{plev|PLEV},{lat|slat-PADDING:nlat+PADDING},{lon|wlon-PADDING:elon+PADDING})
	delete(uwind1)

	vwind1 = fv->va(time|iStart:iLast,plev|:,lat|:,lon|:)
	vwind = vwind1({time|:},{plev|PLEV},{lat|slat-PADDING:nlat+PADDING},{lon|wlon-PADDING:elon+PADDING})
	delete(vwind1)

	delete(uwind_annual)
	delete(vwind_annual)
	delete(uwind_lta)
	delete(vwind_lta)
	;annual
	nyears = nmonths/12
	uwind_annual = new((/nyears,lat2,lon2/), float)
	vwind_annual = new((/nyears,lat2,lon2/), float)
	month_idx = 0
	do i = 0,nyears-1
		uwind_annual(i,:,:) = dim_avg_n_Wrap(uwind(month_idx:month_idx+11,:,:),0)
		vwind_annual(i,:,:) = dim_avg_n_Wrap(vwind(month_idx:month_idx+11,:,:),0)
		month_idx = month_idx + 12
	end do
	;printVarSummary(uwind_annual)
	;printVarSummary(vwind_annual)

	uwind_lta = dim_avg_n_Wrap(uwind_annual,0)
	vwind_lta = dim_avg_n_Wrap(vwind_annual,0)

	;seasonal
	; winter might be tricky because it is across years

	U_MAM = new((/nyears, lat2,lon2/), float)
	V_MAM = new((/nyears, lat2,lon2/), float)
	U_JJA = new((/nyears, lat2,lon2/), float)
	V_JJA = new((/nyears, lat2,lon2/), float)
	U_SON = new((/nyears, lat2,lon2/), float)
	V_SON = new((/nyears, lat2,lon2/), float)

	iYear = 0
	do i = 2, nmonths-1, 12
		; MAM
		U_MAM(iYear,:,:) = dim_avg_n_Wrap(uwind(i:i+2,:,:),0)
		V_MAM(iYear,:,:) = dim_avg_n_Wrap(vwind(i:i+2,:,:),0)
		; JJA
		U_JJA(iYear,:,:) = dim_avg_n_Wrap(uwind(i+3:i+2+3,:,:),0)
		V_JJA(iYear,:,:) = dim_avg_n_Wrap(vwind(i+3:i+2+3,:,:),0)
		; SON
		U_SON(iYear,:,:) = dim_avg_n_Wrap(uwind(i+6:i+2+6,:,:),0)
		V_SON(iYear,:,:) = dim_avg_n_Wrap(vwind(i+6:i+2+6,:,:),0)
		iYear = iYear + 1
	end do
	; average over all the years
	U_MAM_CLM = dim_avg_n_Wrap(U_MAM,0)
	V_MAM_CLM = dim_avg_n_Wrap(V_MAM,0)
	U_JJA_CLM = dim_avg_n_Wrap(U_JJA,0)
	V_JJA_CLM = dim_avg_n_Wrap(V_JJA,0)
	U_SON_CLM = dim_avg_n_Wrap(U_SON,0)
	V_SON_CLM = dim_avg_n_Wrap(V_SON,0)
	

	;test = new(36, float)
	;test(:) = 0
	;test(0:1) = 1
	;do i = 11, 36-1-12, 12
	;	test(i:i+2) = 1
	;end do
	;test(36-2:35) = 1
	;print(test)

	; the DJF is a little different because the last year only has 1 month,
	; so I am going to not include it in the average
	; the first year is going to be only JF
	U_DJF = new((/nyears+1, lat2,lon2/), float)
	V_DJF = new((/nyears+1, lat2,lon2/), float)

	U_DJF(:,:,:) = 0
	V_DJF(:,:,:) = 0

	iYear = 0
	U_DJF(iYear,:,:) = dim_avg_n_Wrap(uwind(0:1,:,:),0)
	V_DJF(iYear,:,:) = dim_avg_n_Wrap(vwind(0:1,:,:),0)
	iYear = 1
	
	do i = 11, nmonths-1-12, 12
		; DJF
		U_DJF(iYear,:,:) = dim_avg_n_Wrap(uwind(i:i+2,:,:),0)
		V_DJF(iYear,:,:) = dim_avg_n_Wrap(vwind(i:i+2,:,:),0)
		iYear = iYear + 1
	end do
	printVarSummary(U_DJF)

	U_DJF_CLM = dim_avg_n_Wrap(U_DJF,0)
	V_DJF_CLM = dim_avg_n_Wrap(V_DJF,0)

	vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag

	plot2 = new(5,graphic)
	
	vcres@tiMainString = "Annual (85000)"
	plot2(0) = gsn_csm_vector_map_ce(wks,uwind_lta,vwind_lta,vcres)
	vcres@tiMainString = "Spring"
	plot2(1) = gsn_csm_vector_map_ce(wks,U_MAM_CLM,V_MAM_CLM,vcres)
	vcres@tiMainString = "Summer"
	plot2(2) = gsn_csm_vector_map_ce(wks,U_JJA_CLM,V_JJA_CLM,vcres)
	vcres@tiMainString = "Autum"
	plot2(3) = gsn_csm_vector_map_ce(wks,U_SON_CLM,V_SON_CLM,vcres)
	vcres@tiMainString = "Winter"
	plot2(4) = gsn_csm_vector_map_ce(wks,U_DJF_CLM,V_DJF_CLM,vcres)

	pidx = 0
	do i = 0,3
		dum(pidx) = gsn_add_polyline(wks, plot1(0), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot1(1), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot1(2), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot1(3), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot1(4), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot2(0), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot2(1), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot2(2), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot2(3), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
		dum(pidx) = gsn_add_polyline(wks, plot2(4), xpts(i:i+1),ypts(i:i+1), respoly)
		pidx = pidx + 1
	end do


	plot = new(10,graphic)
	plot(0) = plot2(0)
	plot(1) = plot1(0)
	plot(2) = plot2(1)
	plot(3) = plot1(1)
	plot(4) = plot2(2)
	plot(5) = plot1(2)
	plot(6) = plot2(3)
	plot(7) = plot1(3)
	plot(8) = plot2(4)
	plot(9) = plot1(4)

	do i=0,9
		draw(plot(i))
		frame(wks)
	end do



	resP = True
	resP@gsnMaximize = True
	gsn_panel(wks,plot,(/5,2/),resP)

	
	exit



;	plot = new(10,graphic)


; res                  = True                     ; plot mods desired
; res@tiYAxisString   = "Moisture Flux (mm/day)"
; res@tiXAxisString   = ""
; res@gsnDraw              = False           ; don't draw
; res@gsnFrame             = False           ; don't advance frame
; res@tmXBMode             = "Explicit"      ; label independently
; res@tmXBValues           = ispan(0,11,1)
; res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
;  "O","N","D"/)
;  res@xyLineThicknessF = 3.0       ; double the thickness of the line
;  res@xyDashPatterns = 0.
;
; time = ispan(0,11,1)
; 
; res@tiMainString     = asciifile+" lat=41N-49N,lon=94W-75W"
; plot(0)  = gsn_csm_xy (wks,time,flux_monthly,res) ; create plot
; delete(res@tmXBMode)
; plot(1)  = gsn_csm_xy (wks,plotyear,flux_annual,res) ; create plot
; asciiwrite(asciidir+asciifile+"_monthly"+asciiext,flux_monthly)
; asciiwrite(asciidir+asciifile+"_annual"+asciiext,flux_annual)
;
;  res@xyLineColors = (/"black","green","blue","red"/)
;
;  res@pmLegendDisplayMode = "Always"
;  res@xyExplicitLegendLabels = (/"south","north","west","east"/)
;  res@pmLegendWidthF = 0.12
;  res@pmLegendHeightF = 0.15
;  res@pmLegendSide = "Top"
;  res@pmLegendOrthogonalPosF = -0.4
;  res@pmLegendParallelPosF = 0.8
;  res@lgPerimOn = False
;
;; plot(2)  = gsn_csm_xy (wks,time,(/flux_monthlyS,flux_monthlyN,flux_monthlyW,flux_monthlyE/),res) ; create plot
;; plot(3)  = gsn_csm_xy (wks,plotyear,(/flux_annualS,flux_annualN,flux_annualW,flux_annualE/),res) ; create plot
;
;   resP                     = True               ; panel options
;  resP@gsnMaximize         = True               ; maximize image
;  resP@gsnPanelLabelBar    = True     
;  
;  gsn_panel(wks,plot,(/2,1/),False)          ; draw: 2-across, 1-down



end
