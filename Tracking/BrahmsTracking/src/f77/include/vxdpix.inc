*  
* $Id: vxdpix.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      REAL VDR(5)
      REAL VXF1(5),VXZF1,VXF2(5),VXZF2,VXF3(5),VXF4(5),VXTR(3)
      REAL VTXD1(3),VTXD2(3),VTXD3(3),VTXD4(3),VTXD5(3),
     & VXTHKF,VXTHKB,VXDPPNT
      COMMON/VXDPIX/VXF1,VXZF1,VXF2,VXZF2,VXF3,VXF4,VXTR,
     &              VTXD1,VTXD2,VTXD3,VTXD4,VTXD5,
     &              VXTHKF,VXTHKB,VXDPPNT, VDR
      SAVE /VXDPIX/

*+KEEP,VXDPIX.
*      REAL VXF1(5),VXZF1,VXF2(5),VXZF2,VXTR(3)
*      REAL VTXD1(3),VTXD2(3),VTXD3(3),VXTHKF,VXTHKB,VXDPPNT
*      COMMON/VXDPIX/VXF1,VXZF1,VXF2,VXZF2,VXTR,
*     &              VTXD1,VTXD2,VTXD3,VXTHKF,VXTHKB,VXDPPNT
*      SAVE /VXDPIX/
