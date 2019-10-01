
def diskmodel(outname='disk', direction='J2000 10h00m00.0s 20d00m00.0s', 
              reffreq='2.8GHz', flux=660000.0, eqradius='16.166arcmin', polradius='16.166arcmin', 
              pangle='21.1deg', index=None, cell='2.0arcsec'):
    ''' Create a blank solar disk model image (or optionally a data cube)
   
        outname       String to use for part of the image and fits file names (default 'disk')
        direction     String specifying the position of the Sun in RA and Dec.  Default
                        means use the standard string "J2000 10h00m00.0s 20d00m00.0s"
        reffreq       The reference frequency to use for the disk model (the frequency at which
                        the flux level applies). Default is '2.8GHz'.
        flux          The flux density, in Jy, for the entire disk. Default is 66 sfu.
        eqradius      The equatorial radius of the disk.  Default is 
                        16 arcmin + 10" (for typical extension of the radio limb)
        polradius     The polar radius of the disk.  Default is 
                        16 arcmin + 10" (for typical extension of the radio limb)
        pangle        The solar P-angle (geographic position of the N-pole of the Sun) in 
                        degrees E of N.  This only matters if eqradius != polradius
        index         The spectral index to use at other frequencies.  Default None means
                        use a constant flux density for all frequencies.
        cell          The cell size (assumed square) to use for the image.  The image size
                        is determined from a standard radius of 960" for the Sun, divided by
                        cell size, increased to nearest power of 512 pixels. The default is '2.0arcsec',
                        which results in an image size of 1024 x 1024.
        Note that the frequency increment used is '325MHz', which is the width of EOVSA bands 
          (not the width of individual science channels)
    '''
    cl.done()
    try:
        aspect = 1.01 #Enlarge the equatorial disk by 1%
        num, unit = eqradius.split('arc')
        diamajor = str(2*float(num)*aspect)+'arc'+unit
        num, unit = polradius.split('arc')
        diaminor = str(2*float(num))+'arc'+unit
        solrad = float(num)
        if unit == 'min':
            solrad *= 60.
    except:
        print 'Radius',radius,'does not have the expected format, number + unit where unit is arcmin or arcsec'
        return
    try:
        num, unit = cell.split('arc')
        cellsize = float(num)
        if unit == 'min':
            cellsize *= 60.
        diskpix = solrad*2/cellsize
        cell_rad = cellsize*pi/180./3600.
    except:
        print 'Cell size',cell,'does not have the expected format, number + unit where unit is arcmin or arcsec'
        return
    # Add 90 degrees to pangle, due to angle definition in addcomponent() -- it puts the majoraxis vertical
    pangle = str(float(pangle.split('deg')[0]) + 90)+'deg'
    mapsize = ((int(diskpix) / 512) + 1)*512
    # Flux density is doubled because it is split between XX and YY
    cl.addcomponent(dir=direction, flux=flux*2, fluxunit='Jy', freq=reffreq, shape='disk', 
                    majoraxis=diamajor, minoraxis=diaminor, positionangle=pangle)
    cl.setrefdirframe(0, 'J2000')
    ia.fromshape(outname+reffreq+'.xim', [mapsize,mapsize,1,1], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad_val = qa.convert(qa.quantity(cell_rad),'rad')['value']
    cs.setincrement([-cell_rad_val,cell_rad_val],'direction')
    epoch, ra, dec = direction.split()
    cs.setreferencevalue([qa.convert(ra,'rad')['value'],qa.convert(dec,'rad')['value']],type="direction")
    cs.setreferencevalue(reffreq,'spectral')
    cs.setincrement('325MHz','spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)
    exportfits(imagename=outname+reffreq+'.xim',fitsimage=outname+reffreq+'.fits',overwrite=True)
    importfits(fitsimage=outname+reffreq+'.fits',imagename=outname+reffreq+'.im',overwrite=True)
    ia.close()
    ia.done()
    cl.close()
    cl.done()
    return
    