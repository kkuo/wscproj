
; save new output variable; Author: Kevin Kuo, Allison Steiner
; ALS modifications 12 Oct 12 based on Zangvil et al., 2004

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 0: Set constants and variables
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)

	; define GL region
	; need to be in float for gc_latlon
	;slat = 40.
	;nlat = 50.
	;elon = 285.
	;wlon = 265.

	slat = 41.
	nlat = 49.
	elon = 285.;-360
	wlon = 266.;-360
	;nmonths = 28*12 ; 1979 to 2007

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12

	asciidir = "."
	asciifile = "Moisture Flux"
	asciiext = ".txt"

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

	asciifnames = (/"/share/kkuo/AR5/scripts/mflux/CCSM4",\
					"/share/kkuo/AR5/scripts/mflux/CNRM-CM5",\
					"/share/kkuo/AR5/scripts/mflux/GFDL-ESM2G", \
					"/share/kkuo/AR5/scripts/mflux/inmcm4",\
					"/share/kkuo/AR5/scripts/mflux/IPSL-CM5A-LR",\
					"/share/kkuo/AR5/scripts/mflux/HadESM2ES",\
					"/share/kkuo/AR5/scripts/mflux/MIROC5",\
					"/share/kkuo/AR5/scripts/mflux/MIROC-ESM",\
					"/share/kkuo/AR5/scripts/mflux/MIROC-ESM-CHEM",\
					"/share/kkuo/AR5/scripts/mflux/MPI-ESM-LR",\
					"/share/kkuo/AR5/scripts/mflux/MRI-CGCM3",\
					"/share/kkuo/AR5/scripts/mflux/NorESM"/)

	husfnames = (/	"/share/kkuo/AR5/CCSM4/keep/hus_Amon_CCSM4_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/CNRM-CM5/hus_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/GFDL-ESM2G/hus.nc", \
					"/share/kkuo/AR5/INMCM4/hus_Amon_inmcm4_1970_2005.nc",\
					"/share/kkuo/AR5/IPSL-CM5A-LR/hus_Amon_IPSL-CM5A-LR_historical_r1i1p1_195001-200512.nc",\
					"/share/kkuo/AR5/HadESM2ES/hus_195912-200511.nc",\
					"/share/kkuo/AR5/MIROC5/hus_Amon_MIROC5_1970_2012.nc",\
					"/share/kkuo/AR5/MIROC-ESM/hus_Amon_MIROC-ESM_1850_2005.nc",\
					"/share/kkuo/AR5/MIROC-ESM-CHEM/hus_Amon_MIROC-ESM-CHEM_1850_2005.nc",\
					"/share/kkuo/AR5/MPI-ESM-LR/hus_1970_to_2005.nc",\
					"/share/kkuo/AR5/MRI-CGCM3/hus_Amon_197001-200512.nc",\
					"/share/kkuo/AR5/NorESM/hus_Amon_NorESM1-M_historical_r1i1p1_195001-200512.nc"/)

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

	nfiles = dimsizes(husfnames)

	rad = 4.0*atan(1.0)/180.0
	re = 6371220.0
	rr = re*rad

	do ifile = 0, nfiles-1
		monFileStrings = asciifnames(ifile)+"monthly.txt"
		annFileStrings = asciifnames(ifile)+"annual.txt"
		;if(isfilepresent(monFileStrings) .and. isfilepresent(annFileStrings))
		;	continue
		;end if

		print(husfnames(ifile))
		cmdstring = (/"ncdump -h ",husfnames(ifile)," | grep units"/)
		cmd = str_concat(cmdstring)
		system(cmd)
		print(husfnames(ifile))
		cmdstring = (/"ncdump -h ",uafnames(ifile)," | grep units"/)
		cmd = str_concat(cmdstring)
		system(cmd)
		print(husfnames(ifile))
		cmdstring = (/"ncdump -h ",vafnames(ifile)," | grep units"/)
		cmd = str_concat(cmdstring)
		system(cmd)

		fq = addfile(husfnames(ifile),"r")
		fu = addfile(uafnames(ifile),"r")
		fv = addfile(vafnames(ifile),"r")

		YYYYMM = cd_calendar(fq->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		;HadESM2-ES data is numbered differently. December appears to be missing from 2005, but it is really not
		if(ismissing(iLast))
			iLast = dimsizes(YYYYMM)-1
		end if

		nq = iLast-iStart+1
		spechum1 = fq->hus(time|iStart:iLast, plev|:, lat|:, lon|:)
		spechum = spechum1({time|:},{plev|:},{lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(spechum1)	
		printMinMax(spechum,True)

		uwind1 = fu->ua(time|iStart:iLast, plev|:, lat|:,lon|:)
		uwind = uwind1({time|:}, {plev|:}, {lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(uwind1)
		printMinMax(uwind,True)

		vwind1 = fv->va(time|iStart:iLast,plev|:,lat|:,lon|:)
		vwind = vwind1({time|:},{plev|:}, {lat|slat-5:nlat+5},{lon|wlon-10:elon+10})
		delete(vwind1)
		printMinMax(vwind,True)

		dims = dimsizes(spechum)
		levs = dims(1)
		lat2 = dims(2)
		lon2 = dims(3)
		;printVarSummary(spechum)
		lats = spechum&lat
		lons = spechum&lon

		; save new output variable
		flux  = new((/nq/),float)
		year  = new((/nq/),integer)
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 3: calc dx and dy for each grid cell 
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; define a 2D dimension for each grid box
		dx = new((/lat2,lon2/),float)
		dy = new((/lat2,lon2/),float)

		dy2 = new(lat2, float)
		dy2(1:lat2-2) = dble2flt(abs(lats(2:lat2-1)-lats(1:lat2-2))*rr)
		dy2(0) = dy2(1)              ;set up dummy values for the first/last row
		dy2(lat2-1) = dy2(lat2-2)

		; calculate dx
		templon = spechum&lon
		;print(templon)
		dlontmp = templon(1)-templon(0)
		print(dlontmp)
		dlon = dlontmp*rr
		; dx is a function of latitude
		do j=0,lon2-1
			do i=1,lat2-2
				dx(i,j) = doubletofloat(dlon*cos(lats(i)*rad))
			end do
		end do
		dx(0,:) = dx(1,:)           ; set up dummy vars for first and last row
		dx(lat2-1,:) = dx(lat2-2,:)

		dx!0 = "lat"
		dx!1 = "lon"
		dx&lat = lats
		dx&lon = lons
		dy!0 = "lat"
		dy!1 = "lon"
		dy&lat = lats
		dy&lon = lons
		; conform dy to the same grid (dy does not change with dx)
		dy = conform_dims(dimsizes(dx),dy2,0)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 4: calculate dp for each pressure level
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		;fl = addfile("shum.mon.mean.nc","r")
		;pres = fl->level
		pres       = spechum&plev
		plev       = pres
		plev@units = "Pa"
		ptop       = min(plev)
		psfc       = 1013.*100.           ; use dummy value until get the 2D fields
		psfc@units = "Pa"
		dp         = dpres_plevel(plev,psfc,ptop,0)      ; calculate pressure thicknesses
		dp_3d      = conform_dims(dimsizes(uwind),dp,1)  ; need to make the dp array into 3D (use conform)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 5:  Moisture Flux Divergence Calculation
			;   Use Green's theorem as in Zangvil et al., 2004, Eqn 2
		;   Inner integral(dl) - q*v*dl [kg H2O/kg air * m/s * m = kg H2O m2/(kg air s)] 
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

		; conform dims for the dx and dy to make it 3D 
		newdx = conform(spechum,dx,(/2,3/))
		newdy = conform(spechum,dy,(/2,3/))
		copy_VarMeta(spechum,newdx)
		copy_VarMeta(spechum,newdy)
		;printVarSummary(newdx)

		;print(newdx(0,{slat},{wlon:elon}))

		; copy metadata
		quflx = spechum								
		qvflx = spechum 
		; uwind moisture flux divergence (kg m2)/(kg s) at each pressure level
		; for uwind, use the dy length scale
		quflx = spechum*uwind*newdy				
		; vwind moisture flux divergence (kg m2)/(kg s) at each pressure level
		; for vwind, use the dx length scale
		qvflx = spechum*vwind*newdx				

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		;   Outer integral(dp) - inner integral*dp   [kg H2O m2/(kg air s) * kg/m/s2 = kg H2O m/s3]
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

		; new units will be kg m3/s3
		quflx_dp = quflx   ; copy metadata
		qvflx_dp = qvflx   ; copy metadata

		quflx_dp = quflx*dp_3d
		qvflx_dp = qvflx*dp_3d

		; integrate vertically by summing over all heights
		quflx_int = dim_sum_n_Wrap(quflx_dp,1)
		qvflx_int = dim_sum_n_Wrap(qvflx_dp,1) 

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;  Now sum over each boundary (E,W,N,S) = kg H2O m/s3
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ;print("Wbdy (kg/m3):"+quflx_int({slat:nlat},{wlon}))
        ;print("Ebdy (kg/m3):"+quflx_int({slat:nlat},{elon}))
        ;print("Sbdy (kg/m3):"+quflx_int({slat},{wlon:elon}))
        ;print("Nbdy (kg/m3):"+quflx_int({nlat},{wlon:elon}))
        

        W = dim_sum_Wrap(quflx_int(:,{slat:nlat},{wlon}))
        E = dim_sum_Wrap(quflx_int(:,{slat:nlat},{elon}))
        N = dim_sum_Wrap(qvflx_int(:,{nlat},{wlon:elon}))
        S = dim_sum_Wrap(qvflx_int(:,{slat},{wlon:elon}))

        ;; u+ in the west to east direction (+u on E means out flux)
		;; v+ in the south to north direction (+v on N means out flux)
        allsides = (E-W+N-S)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;  Final integral step: divide by Ag = [kg H2O m/s3 /m2 / (m/s2) = kg H2O/m2/s]
		;  Convert units to mm/day (above/rhoL*1000 mm/m) = 
        ;               [kg H2O/m2/s * (m3/kg H2O) * 1000 mm/m *86400 s/day]
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        area_grid = dx   ; copy meta?
        area_grid = dx*dy
        ; sum area of all individual grid boxes
        xyarea = dim_sum_Wrap(area_grid(:,{wlon:elon}))
        xyarea2 = dim_sum_Wrap(xyarea({nlat:slat}))
		xyarea3 = 1327828198707.01 ;m from matlab's areaquad function

        ;print(allsides)
        print(xyarea2)

        flux = allsides/xyarea2/g/rhowat*1000.*86400. 
		;printVarSummary(flux)


		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
		; Step 8 find the montly climatology
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

		; climatology of all years
		flux_monthly = new(12,float)
		do i=0,11
		 flux_monthly(i) = dim_avg_Wrap(flux(i:nq-1:12))
		end do
		   
		; annual value (for comparison wiht Li et al 2010, Fig 6)
		flux_annual = new(nq/12,float)
		count = 0
		do i=0,nq/12-1
		  flux_annual(i) = dim_avg_Wrap(flux(count:count+11))
		 count = count + 12
		end do

		;plotyear = ispan(1979,1979-1+nmonths/12,1)
		;plotyear = new(nq/12,integer)
		;do i=0,nq-1,12
		;   plotyear(i/12) = year(i) 
		;end do

		asciiwrite(monFileStrings,flux_monthly)
		asciiwrite(annFileStrings,flux_annual)

		delete(YYYYMM)
		delete(iStart)
		delete(iLast)
		delete(flux_monthly)
		delete(flux_annual)
		delete(spechum)
		delete(uwind)
		delete(vwind)
		delete(dims)
		delete(levs)
		delete(lat2)
		delete(lon2)
		delete(lats)
		delete(lons)
		delete(flux)
		delete(year)
		delete(dx)
		delete(dy)
		delete(dy2)
		delete(templon)
		delete(pres)
		delete(plev)
		delete(ptop)
		delete(psfc)
		delete(dp)
		delete(dp_3d)
		delete(newdx)
		delete(newdy)
		delete(quflx)
		delete(qvflx)
		delete(quflx_dp)
		delete(qvflx_dp)
		delete(quflx_int)
		delete(qvflx_int)
        delete(W)
        delete(E)
        delete(N)
        delete(S)
        delete(allsides)
        delete(area_grid)
        delete(xyarea)


	end do

	

		





;************************************************
; plotting parameters
;************************************************
 ;wks   = gsn_open_wks ("pdf","narr_mflux_regrid")                ; open workstation

 ;plot = new(4,graphic)



 ;res                  = True                     ; plot mods desired
 ;res@tiYAxisString   = "Moisture Flux (mm/day)"
 ;res@tiXAxisString   = ""
 ;res@gsnDraw              = False           ; don't draw
 ;res@gsnFrame             = False           ; don't advance frame
 ;res@tmXBMode             = "Explicit"      ; label independently
 ;res@tmXBValues           = ispan(0,11,1)
 ;res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S",\
 ; "O","N","D"/)
 ; res@xyLineThicknessF = 3.0       ; double the thickness of the line
 ; res@xyDashPatterns = 0.

 ;time = ispan(0,11,1)
 ;
 ;res@tiMainString     = "NARR regrid lat=41N-49N,lon=94W-75W"
 ;plot(0)  = gsn_csm_xy (wks,time,flux_monthly,res) ; create plot
 ;delete(res@tmXBMode)
 ;plot(1)  = gsn_csm_xy (wks,plotyear,flux_annual,res) ; create plot

  ;res@xyLineColors = (/"black","green","blue","red"/)

  ;res@pmLegendDisplayMode = "Always"
  ;res@xyExplicitLegendLabels = (/"south","north","west","east"/)
  ;res@pmLegendWidthF = 0.12
  ;res@pmLegendHeightF = 0.15
  ;res@pmLegendSide = "Top"
  ;res@pmLegendOrthogonalPosF = -0.4
  ;res@pmLegendParallelPosF = 0.8
  ;res@lgPerimOn = False

; plot(2)  = gsn_csm_xy (wks,time,(/flux_monthlyS,flux_monthlyN,flux_monthlyW,flux_monthlyE/),res) ; create plot
; plot(3)  = gsn_csm_xy (wks,plotyear,(/flux_annualS,flux_annualN,flux_annualW,flux_annualE/),res) ; create plot

	;resP                     = True               ; panel options
	;resP@gsnMaximize         = True               ; maximize image
	;resP@gsnPanelLabelBar    = True     
  
  ;gsn_panel(wks,plot,(/2,1/),False)          ; draw: 2-across, 1-down



end
