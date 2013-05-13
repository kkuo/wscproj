load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	g = 9.8        ; gravitational constant (m/s2)
	rhowat = 1000. ; density of water (kg/m3)

	; define GL region
	; need to be in float for gc_latlon
	slat = 40.
	nlat = 50.
	elon = 285.
	wlon = 265.

	;these are what Li uses
	;slat = 41.
	;nlat = 49.
	;elon = 285.;-360
	;wlon = 266.;-360

	yrStart = 1979
	yrLast = 2005
	nq = (yrLast-yrStart+1)*12

	asciidir = "."
	asciifile = "EvapCont"
	asciiext = ".txt"


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	; Step 1:  Read in data and get dimensions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

	fnames = (/	"/share/kkuo/AR5/CCSM4/keep/evspsbl_Amon_CCSM4_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/CNRM-CM5/evspsbl_Amon_CNRM-CM5_historical_r1i1p1_195001-200512.nc",\
				"/share/kkuo/AR5/GFDL-ESM2G/evspsbl_196101-200512.nc", \
				"/share/kkuo/AR5/INMCM4/evspsbl_Amon_inmcm4_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/IPSL-CM5A-LR/evspsbl_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MIROC5/evspsbl_Amon_MIROC5_historical_r1i1p1_185001-201212.nc",\
				"/share/kkuo/AR5/MIROC-ESM/evspsbl_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MIROC-ESM-CHEM/evspsbl_Amon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MPI-ESM-LR/evspsbl_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/MRI-CGCM3/evspsbl_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512.nc",\
				"/share/kkuo/AR5/HadESM2ES/evspsbl_195912-200511.nc",\
				"/share/kkuo/AR5/NorESM/evspsbl_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc"/)
	shortnames = (/	"CCSM4",\
				"CNRM-CM5",\
				"GFDL-ESM2G", \
				"INMCM4",\
				"IPSL-CM5A-LR",\
				"MIROC5",\
				"MIROC-ESM",\
				"MIROC-ESM-CHEM",\
				"MPI-ESM-LR",\
				"MRI-CGCM3",\
				"HadGEM2-ES",\
				"NorESM"/)

	nfiles = dimsizes(fnames)
	wks   = gsn_open_wks ("pdf",asciifile)                ; open workstation
	gsn_define_colormap(wks,"rainbow")

	plot = new(nfiles*2,graphic) ; 2 - one for JJA and the other DJF
	res                  = True                     ; plot mods desired
	res@gsnDraw              = False           ; don't draw
	res@gsnFrame             = False           ; don't advance frame
	res@gsnMaximize = True
	res@cnFillOn = True
	res@tiMainOn = True
	res@lbLabelBarOn = False
	res@cnLinesOn = False
	res@cnLineLabelsOn = True
	res@cnLevelSelectionMode = "ManualLevels"
	res@cnLevelSpacingF = .25
	res@cnMinLevelValF = 0.
	res@cnMaxLevelValF = 5.
	res@gsnSpreadColors = True
	res@gsnAddCyclic = False
	res@mpMinLatF = slat - 1
	res@mpMaxLatF = nlat + 1
	res@mpMinLonF = wlon - 1
	res@mpMaxLonF = elon + 1


	

	do ifile = 0, nfiles-1
		print(fnames(ifile))
		cmdstring = (/"ncdump -h ",fnames(ifile)," | grep pr:units"/)
		cmd = str_concat(cmdstring)
		system(cmd)

		isHadESM = False

		fp = addfile(fnames(ifile),"r")
		YYYYMM = cd_calendar(fp->time,-1)
		iStart = ind(YYYYMM.eq.(yrStart*100+1))
		iLast = ind(YYYYMM.eq.(yrLast*100+12))
		if(ismissing(iLast))
			isHadESM = True
			iLast = dimsizes(YYYYMM)-12
		end if

		pr1 = fp->evspsbl(time|iStart:iLast, lat|:, lon|:)
		pr = pr1({time|:},{lat|slat:nlat},{lon|wlon:elon})
		dims_pr = dimsizes(pr)
		nmonths = dims_pr(0)
		nyears = nmonths/12
		

		iyear = 0
		prjja = new((/nyears,dims_pr(1),dims_pr(2)/),float)
		;print(nmonths)
		;print(nyears)
		do imon=5,nmonths-1,12
			;print(imon)
			prjja(iyear,:,:) = dim_avg_n_Wrap(pr(imon:imon+2,:,:),0)
			iyear = iyear+1
		end do
		prjja2 = dim_avg_n(prjja,0)*86400
		prjja2!0 = "lat"
		prjja2!1 = "lon"
		prjja2&lat = pr&lat
		prjja2&lon = pr&lon
		prjja2&lat@units = "degrees_north"
		prjja2&lon@units = "degrees_east"
		printMinMax(prjja2,True)

		res@tiMainString = shortnames(ifile)
		plot(ifile*2) = gsn_csm_contour_map(wks,prjja2,res)

		;DJF
		prdjf = new((/nyears, dims_pr(1),dims_pr(2)/), float)

		iYear = 0
		prdjf(iYear,:,:) = dim_avg_n_Wrap(pr(0:1,:,:),0)
		iYear = 1
		
		do i = 11, nmonths-1-12, 12
			prdjf(iYear,:,:) = dim_avg_n_Wrap(pr(i:i+2,:,:),0)
			iYear = iYear + 1
		end do
		prdjf2 = dim_avg_n(prdjf,0)*86400
		prdjf2!0 = "lat"
		prdjf2!1 = "lon"
		prdjf2&lat = pr&lat
		prdjf2&lon = pr&lon
		prdjf2&lat@units = "degrees_north"
		prdjf2&lon@units = "degrees_east"
		printMinMax(prdjf2,True)
		plot(ifile*2+1) = gsn_csm_contour_map(wks,prdjf2,res)


		delete(YYYYMM)
		delete(iStart)
		delete(iLast)
		delete(pr1)
		delete(pr)
		delete(prjja)
		delete(prjja2)
		delete(prdjf)
		delete(prdjf2)
	end do
	resP = True
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True
	;res@lbAutoManage = True
	resP@lbLabelFontHeightF = 0.01
	resP@lbLabelStride = 2
	nfiles2 = nfiles/2
	offset = mod(nfiles,2)
	print(nfiles)
	print(nfiles2)
	gsn_panel(wks,plot(0:nfiles2*2-1),(/nfiles2,2/),resP)
	gsn_panel(wks,plot(nfiles2*2:),(/nfiles2+offset,2/),resP)

;************************************************
; plotting parameters
;************************************************





end
