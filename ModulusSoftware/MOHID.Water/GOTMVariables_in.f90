!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License, 
!as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

      type(T_GOTM), pointer     :: ObjGOTM

      ! Variables in ObjGOTM
      double precision,  dimension(:), pointer    :: cmue1,cmue2
      double precision,  dimension(:), pointer    :: tkeo
      double precision,  dimension(:), pointer    :: as,an
      double precision,  dimension(:), pointer    :: xRf 

      ! Variables in ObjGOTM%Export
      double precision,  dimension(:), pointer    :: tke,eps,L
      double precision,  dimension(:), pointer    :: num,nuh

      ! Variables in ObjGOTM%Trid
      double precision,  dimension(:), pointer    :: au,bu,cu,du,ru,qu

      ! Variables in ObjGOTM%Parameters
      double precision  const_num,const_nuh,k_min,L_min,eps_min
      double precision kappa,Prandtl0,cm0,cm_craig,cw,galp
      double precision ce1,ce2,ce3minus,ce3plus,sig_k
      double precision sl,e1,e2,e3
      double precision a1,a2,b1,b2,c2,c3,qeghmax,qeghmin,qeghcrit   
      double precision alpha,klimiw,rich_cr,numiw,nuhiw,numshear
      double precision c1,cde,craig_m,sig_e0,sig_e1
      integer turb_method,tke_method,len_scale_method,stab_method,MY_length
      integer iw_model   
      logical craig_banner,length_lim
      logical qesmooth,flux_bdy
      
      ! Variables in ObjGOTM%Export
      tke   =>   ObjGOTM%Export%tke 
      eps   =>   ObjGOTM%Export%eps 
      L     =>   ObjGOTM%Export%L   
      num   =>   ObjGOTM%Export%num 
      nuh   =>   ObjGOTM%Export%nuh 

      !Variables in ObjGOTM%Trid
      au    =>   ObjGOTM%Trid%au       
      bu    =>   ObjGOTM%Trid%bu       
      cu    =>   ObjGOTM%Trid%cu     
      du    =>   ObjGOTM%Trid%du 
      ru    =>   ObjGOTM%Trid%ru     
      qu    =>   ObjGOTM%Trid%qu        

      ! Variables in ObjGOTM
      cmue1        => ObjGOTM%cmue1
      cmue2        => ObjGOTM%cmue2
      tkeo         => ObjGOTM%tkeo
      as           => ObjGOTM%as
      an           => ObjGOTM%an
      xRf          => ObjGOTM%xRf
      
      ! Variables in ObjGOTM%Parameters
      c1           = ObjGOTM%Parameters%c1
      cde          = ObjGOTM%Parameters%cde
      craig_m      = ObjGOTM%Parameters%craig_m
      sig_e0       = ObjGOTM%Parameters%sig_e0
      sig_e1       = ObjGOTM%Parameters%sig_e1

      turb_method  = ObjGOTM%Parameters%turb_method
      tke_method   = ObjGOTM%Parameters%tke_method
      len_scale_method    = ObjGOTM%Parameters%len_scale_method
      stab_method  = ObjGOTM%Parameters%stab_method
      craig_banner = ObjGOTM%Parameters%craig_banner
      length_lim   = ObjGOTM%Parameters%length_lim
      const_num    = ObjGOTM%Parameters%const_num
      const_nuh    = ObjGOTM%Parameters%const_nuh
      k_min        = ObjGOTM%Parameters%k_min
      L_min        = ObjGOTM%Parameters%L_min
      eps_min      = ObjGOTM%Parameters%eps_min
      kappa        = ObjGOTM%Parameters%kappa
      Prandtl0     = ObjGOTM%Parameters%Prandtl0
      cm0          = ObjGOTM%Parameters%cm0
      cm_craig     = ObjGOTM%Parameters%cm_craig
      cw           = ObjGOTM%Parameters%cw
      galp         = ObjGOTM%Parameters%galp
      ce1          = ObjGOTM%Parameters%ce1
      ce2          = ObjGOTM%Parameters%ce2
      ce3minus     = ObjGOTM%Parameters%ce3minus
      ce3plus      = ObjGOTM%Parameters%ce3plus
      sig_k        = ObjGOTM%Parameters%sig_k
      flux_bdy     = ObjGOTM%Parameters%flux_bdy
      sl           = ObjGOTM%Parameters%sl
      e1           = ObjGOTM%Parameters%e1
      e2           = ObjGOTM%Parameters%e2
      e3           = ObjGOTM%Parameters%e3
      MY_length    = ObjGOTM%Parameters%MY_length
      a1           = ObjGOTM%Parameters%a1
      a2           = ObjGOTM%Parameters%a2
      b1           = ObjGOTM%Parameters%b1
      b2           = ObjGOTM%Parameters%b2
      c2           = ObjGOTM%Parameters%c2
      c3           = ObjGOTM%Parameters%c3
      qesmooth     = ObjGOTM%Parameters%qesmooth
      qeghmax      = ObjGOTM%Parameters%qeghmax
      qeghmin      = ObjGOTM%Parameters%qeghmin
      qeghcrit     = ObjGOTM%Parameters%qeghcrit
      iw_model     = ObjGOTM%Parameters%iw_model
      alpha        = ObjGOTM%Parameters%alpha
      klimiw       = ObjGOTM%Parameters%klimiw
      rich_cr      = ObjGOTM%Parameters%rich_cr
      numiw        = ObjGOTM%Parameters%numiw
      nuhiw        = ObjGOTM%Parameters%nuhiw
      numshear     = ObjGOTM%Parameters%numshear

!----------------------------------------------------------------------------------------------------------
!Copyright (C) 2000 - GOTM code.
!----------------------------------------------------------------------------------------------------------
