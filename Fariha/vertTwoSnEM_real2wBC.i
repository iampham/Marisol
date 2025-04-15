
# TODO: Define kappaSnGrain1, kappaSnCu, kappaSnSn, kappaGrain1Cu

[Mesh]
  type = FileMesh
  file = two_vert_Cu_noIMC.msh
[]

[Variables]

  # concentration

  [./psi] # electric potential
  order = FIRST
  family = LAGRANGE
  #scaling = 1e3
[../]  

  [c]
    order = FIRST
    family = LAGRANGE
  []

  # order parameter 1
  [etaSn1]
    order = FIRST
    family = LAGRANGE
  []

  # order parameter 2
  [etaGrain1]
    order = FIRST
    family = LAGRANGE
  []

  # order parameter 3
  [etaSn2]
    order = FIRST
    family = LAGRANGE
  []

  [etaCu]
    order = FIRST
    family = LAGRANGE
  []

  # phase concentration 1
  [cSn1]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e-3
  []

  # phase concentration 2
  [cGrain1]
    order = FIRST
    family = LAGRANGE
  []

  # phase concentration 3
  [cSn2]
    order = FIRST
    family = LAGRANGE
  []

  [cCu]
    order = FIRST
    family = LAGRANGE
  []


  # Lagrange multiplier
  [lambda]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []


[]

[AuxVariables]

  [detaSn1dt]
  []

  [acbulkf]
  []

  [acbulkc]
  []

  [acinterface]
  []

  [lambdakernel]
  []

  [dcdt]
  []

[]





[BCs]

  [./Top_psi_neumann]
      type = NeumannBC
      variable = 'psi'
      boundary = 'Top'
      value = 5e-5
  [../]
  
  
  [./Bottom_psi_neumann]
      type = NeumannBC
      variable = 'psi'
      boundary = 'Bottom'
      value = -5e-5
  [../]


[]

[ICs]

  [psi]
    variable = psi
    type = ConstantIC
    value = 0
  []

  [etaSn1]
    variable = etaSn1
    type = ConstantIC
    value = 1.0
    block = 'Sn1' 
  []
  [etaGrain1]
    variable = etaGrain1
    type = ConstantIC
    value = 1.0
    block = 'IMC'
  []
  [etaSn2]
    variable = etaSn2
    type = ConstantIC
    value = 1.0
    block = 'Sn2'
  []

  [etaCu] 
    variable = etaCu
    type = ConstantIC
    value = 1
    block = 'Cu'
  []

  [c_Sn]
    variable = c
    type = ConstantIC
    value = 0.05
    block = 'Sn1 Sn2'
  []
  [c_Cu6Sn5]
    variable = c
    type = ConstantIC
    value = 0.6
    block = 'IMC'
    # block = 'Cu6Sn5_1'
  []

  [c_Cu]
    variable = c
    type = ConstantIC
    value = 0.9
    block = 'Cu'
  []


[]

[Materials]

  [initialTemperature]
    type = GenericConstantMaterial
    prop_names = 'T_initial'
    prop_values = '293'
  []

  [molarVolume] 
    type = GenericConstantMaterial
    prop_names = 'Vm          scale    scaleEM' 
    prop_values = '1.629e-5   5e-6     0'
                  #m^3/mol    #converts to pJ/um^3... factor of 1e-6 multiplied by  5 for arbitrary scaling
  []
  # simple toy free energies
  [fSn1] # J/m^3 scaled by factor at end to become pJ/um^3
    type = DerivativeParsedMaterial
    f_name = FSn1
    args = 'cSn1'
    material_property_names = 'Vm scale'
    # function = '2e2*(cSn-0.9)^2+50'
    #multiply the below function by 1e4 to get true result
    # function = '(4.206e2/2*((1-cSn)-0.99)^2+7.168e-1*((1-cSn)-0.99)-1.53)*2e2'
    # function = '(1.9565e0*(1-cSn)^2 -3.3694e0*(1-cSn) - 1.3980e0)*2e2'
    # function = '(0.5641e1*(1-cSn)^4 -1.4942e1*(1-cSn)^3 + 1.4706e1*(1-cSn)^2-0.7204*(1-cSn)-0.1125)*2e1'
    #function = '(0.5641e5*(1-cSn1)^4 -1.4942e5*(1-cSn1)^3 + 1.4706e5*(1-cSn1)^2-0.7204e5*(1-cSn1)-0.1125e5)/Vm*scale'
    function = '(0.5513e5*(1-cSn1)^4 -1.4666e5*(1-cSn1)^3 + 1.44426e5*(1-cSn1)^2-0.7061e5*(1-cSn1)-0.1123e5)/Vm*scale'
    derivative_order = 2
  []
   [fSn2] # J/m^3 scaled by factor at end to become pJ/um^3
    type = DerivativeParsedMaterial
    f_name = FSn2
    args = 'cSn2'
    material_property_names = 'Vm scale'
    # function = '2e2*(cSn-0.9)^2+50'
    #multiply the below function by 1e4 to get true result
    # function = '(4.206e2/2*((1-cSn)-0.99)^2+7.168e-1*((1-cSn)-0.99)-1.53)*2e2'
    # function = '(1.9565e0*(1-cSn)^2 -3.3694e0*(1-cSn) - 1.3980e0)*2e2'
    # function = '(0.5641e1*(1-cSn)^4 -1.4942e1*(1-cSn)^3 + 1.4706e1*(1-cSn)^2-0.7204*(1-cSn)-0.1125)*2e1'
    #function = '(0.5641e5*(1-cSn2)^4 -1.4942e5*(1-cSn2)^3 + 1.4706e5*(1-cSn2)^2-0.7204e5*(1-cSn2)-0.1125e5)/Vm*scale'
    function = '(0.5513e5*(1-cSn2)^4 -1.4666e5*(1-cSn2)^3 + 1.44426e5*(1-cSn2)^2-0.7061e5*(1-cSn2)-0.1123e5)/Vm*scale'
    derivative_order = 2
  []
  
  [fGrain1] # J/m^3 scaled by factor at end 
    type = DerivativeParsedMaterial
    f_name = FGrain1
    args = 'cGrain1'
    material_property_names = 'Vm scale'   
    # function = '2e2*(cGrain1-0.4)^2-100'
    #multiply the below function by 1e4 to get true result
    # function = '(2.0e1/2*((1-cGrain1)-0.41753)^2-6.9892e-1*((1-cGrain1)-0.41753)-1.92)*2e2'
    # function = '(2.0e1*(1-cGrain1)^2-1.74e1*(1-cGrain1)+0.0762e1)*2e1'
    # function = '(2.0e5*(1-cGrain1)^2-1.74e5*(1-cGrain1)+0.0762e5)/Vm*scale'
    function =  '(8.0e5*(1-cGrain1)^2-6.96e5*(1-cGrain1)+1.2117e5)/Vm*scale'
    derivative_order = 2
  []
  
  [fCu]# J/m^3 scaled by factor at end
    type = DerivativeParsedMaterial
    f_name = FCu
    args = 'cCu'
    material_property_names = 'Vm scale'
    # function = '2e2*(cCu-0.1)^2'  
    #multiply the below function by 1e4 to get true result
    # function = '(1.0133e1/2*((1-cCu)-0.105)^2-2.11*((1-cCu)-0.105)-1.28)*2e2'
    # function = '(2.215e0*(1-cCu)^2-2.801e0*(1-cCu)-2.0661e0)*2e2'
    # function = '(1.4397*(1-cCu)^4 -5.2174*(1-cCu)^3 + 7.5859*(1-cCu)^2-4.6290*(1-cCu)-1.9238)*2e1'
    #function = '(1.4397e4*(1-cCu)^4 -5.2174e4*(1-cCu)^3 + 7.5859e4*(1-cCu)^2-4.6290e4*(1-cCu)-1.9238e4)/Vm*scale'
    function = '(1.3021e4*(1-cCu)^4 -4.9421e4*(1-cCu)^3 + 7.3054e4*(1-cCu)^2-4.4861e4*(1-cCu)-1.9222e4)/Vm*scale'
    derivative_order = 2

  []

  

  [hSn1]
    type = DerivativeParsedMaterial
    f_name = hSn1
    function = 'etaSn1^2/(etaSn1^2+etaGrain1^2+etaSn2^2+etaCu^2)'
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
  []

  [hGrain1]
    type = DerivativeParsedMaterial
    f_name = hGrain1
    function = 'etaGrain1^2/(etaSn1^2+etaGrain1^2+etaSn2^2+etaCu^2)'
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
  []

  [hSn2]
    type = DerivativeParsedMaterial
    f_name = hSn2
    function = 'etaSn2^2/(etaSn1^2+etaGrain1^2+etaSn2^2+etaCu^2)'
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
  []


  [hCu]
    type = DerivativeParsedMaterial
    f_name = hCu
    function = 'etaCu^2/(etaSn1^2+etaGrain1^2+etaSn2^2+etaCu^2)'
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
    
  []
  [hintSn]
    type = DerivativeParsedMaterial
    f_name = hintSn
    function = 'hSn1*(hGrain1) + hSn2*hGrain1'
    args = 'etaSn1 etaGrain1 etaSn2'
    material_property_names = 'hSn1 hGrain1 hSn2'
    
  []

  [hintCu]
    type = DerivativeParsedMaterial
    f_name = hintCu
    function = 'hCu*(hGrain1)'
    args = 'etaGrain1 etaCu'
    material_property_names = 'hCu hGrain1 hSn2'
    
  []

  [hint1]
    type = DerivativeParsedMaterial
    f_name = hint1
    function = 'hSn1*hCu'
    args = 'etaSn1 etaSn2 etaCu'
    material_property_names = 'hCu hSn1 hSn2'
    outputs = exodus
    
  []


  [hint2]
    type = DerivativeParsedMaterial
    f_name = hint2
    function = 'hSn2*hCu'
    args = 'etaSn1 etaSn2 etaCu'
    material_property_names = 'hCu hSn2'
    outputs = exodus
    
  []

   [hGB]
     type = DerivativeParsedMaterial
     f_name = hGB
     function = 'hSn1*hSn2'
     args = 'etaSn1 etaSn2'
     material_property_names = 'hSn1 hSn2'
     outputs = exodus
    
   []

  
 

   # Barrier functions for each phase
  [gSn1] 
    type = DerivativeParsedMaterial
    f_name = gSn1
    function = '1.6*(etaSn1^2*etaSn2^2) + 1.25*(etaSn1^2*etaGrain1^2) + 4.8*(etaSn1^2 * etaCu^2)'
    #function = 'w*abs(etaSn1*etaSn2) + w*abs(etaSn1*etaGrain1) + w*abs(etaSn1*etaCu)'
    #function = 'w*hSn1*hSn2 + w*hSn1*hGrain1 + w*hSn1*hCu'
    args = 'etaSn1 etaSn2 etaGrain1 etaCu'
    outputs = exodus
  []
  [gSn2] 
    type = DerivativeParsedMaterial
    f_name = gSn2
    function = '1.6*(etaSn1^2 * etaSn2^2) + 1.25*(etaSn2^2 * etaGrain1^2) + 4.8*(etaSn2^2 * etaCu^2)'
    args = 'etaSn1 etaSn2 etaGrain1 etaCu'
    # material_property_names = 'wSnSn wSnGrain1 wSnCu'
    #function = 'w*hSn1*hSn2 + w*hSn2*hGrain1 + w*hSn2*hCu'
    #material_property_names = 'w hSn1 hSn2 hGrain1 hCu'
  []
  [gGrain1]
    type = DerivativeParsedMaterial
    f_name = gGrain1
    function = '1.25*(etaGrain1^2 * etaSn1^2) + 1.25*(etaGrain1^2 * etaSn2^2) + 5*(etaGrain1^2 * etaCu^2)'
    #function = 'w*abs(etaGrain1*etaSn1) + w*abs(etaGrain1*etaSn2) + w*abs(etaGrain1*etaCu) '
    args = 'etaSn1 etaSn2 etaGrain1 etaCu'
  []

  [gCu] 
    type = DerivativeParsedMaterial
    f_name = gCu
    function = '5*(etaCu^2 * etaGrain1^2) + 4.8*(etaSn1^2 * etaCu^2) + 4.8*(etaSn2^2 *etaCu^2)'
    args = 'etaSn1 etaSn2 etaGrain1 etaCu'
    # material_property_names = 'wSnCu wGrain1Cu'
    #function = 'w*hCu*hSn1 + w*hCu*hGrain1 + w*hSn2*hCu'
    #material_property_names = 'w hSn1 hSn2 hGrain1 hCu'
  []
  
  
  
  [kappa]
    type = ParsedMaterial
    f_name = kappa
    function = '0.04*(hSn1*hSn2 * 0.16 + (hSn1+hSn2)*hGrain1*0.125 + hCu*hGrain1*0.5 + (hSn1+hSn2)*hCu*0.48)*5e0'
    # function = '(hGrain1*hGrain2 * 1.65 + hSn*(hGrain1+hGrain2)*0.66 + hCu*(hGrain1+hGrain2)*1.65)*1e1'
    # function = 1
    material_property_names = 'hSn1 hGrain1 hSn2 hCu'
  []


  [num_phases]
    type = ParsedMaterial
    f_name = num_phases
    function = '2*(if(hCu>0.1,1,0) + if(hSn1>0.1,1,0) + if(hGrain1>0.1,1,0) + if(hSn2>0.1,1,0))'
    material_property_names = 'hCu hSn1 hGrain1 hSn2'
    outputs = exodus
  []
  
  # Note that the cu6sn5/cu reaction must be faster than the GB reaction. If it isn't we get weird tails.
  [L]
    type = ParsedMaterial
    f_name = L
    function = '(hSn1*hSn2 * 0.5 + (hSn1+hSn2)*(hGrain1)*1 + hCu*(hGrain1)*0.5)*100'
    # function = 1
    material_property_names = 'hSn1 hGrain1 hSn2 hCu num_phases'
    args = 'etaCu etaSn1 etaSn2'
    outputs = exodus
  []


  ####################################
  # Diffusivity tensors (using the eigenstrain class which is basically weighting them in MOOSE)
  ###################################

  # Units of um^2/s, assume c-axis for thermal aging.
  [DSn1] # C oriented 90 degrees with y-axis = 
    type = ConstantAnisotropicMobility
    tensor = '5.37e-1 0 0 0 1.42e1 0 0 0 0'
    M_name = DSn1
  []
  [DSn2] # C oriented 90 degrees with y-axis = 
    type = ConstantAnisotropicMobility
    tensor = '1.42e1 0 0 0 5.37e-1 0 0 0 0'
    M_name = DSn2
  []

  [DCu6Sn5]
    type = ConstantAnisotropicMobility
    tensor = '14.2e-4 0 0 0 14.2e-4 0 0 0 0'
    M_name = DCu6Sn5
  []

  [DCu]
    type = ConstantAnisotropicMobility
    tensor = '2.82e-15 0 0 0 2.82e-15 0 0 0 0'
    M_name = DCu
  []

  [DGB]
     type = ConstantAnisotropicMobility
     tensor = '1.42e2 0 0 0 1.42e2 0 0 0 0'
     M_name = DGB
   []

  [DintSn]
    type = ConstantAnisotropicMobility
    tensor = '0.71e1 0 0 0 0.71e1 0 0 0 0'
    M_name = DintSn
  []

  # used to be 0.005
  [DintCu]
    type = ConstantAnisotropicMobility
    tensor = '0.71e1 0 0 0 0.71e1 0 0 0 0'
    M_name = DintCu
  []

  [Dint1]
    type = ConstantAnisotropicMobility
    tensor = '1.58e0 0 0 0 7.64e0 0 0 0 0'
    # tensor = '0 0 0 0 0 0 0 0 0'
    M_name = Dint1
  []

  [Dint2]
    type = ConstantAnisotropicMobility
    tensor = '7.64e0 0 0 0 1.58e0 0 0 0 0'
    # tensor = '0 0 0 0 0 0 0 0 0'
    M_name = Dint2
  []


  [Adjustment]
    type = ParsedMaterial
    f_name = adjustment
    function = '50'
  []

  [a]
    type = ParsedMaterial
    f_name = a
    function = 'if(hSn1>=0.8, 1, 0)'
    material_property_names = 'hSn1'
  []

  [ahSn1]
    type = ParsedMaterial
    f_name = ahSn1
    function = a*hSn1
    material_property_names = 'a hSn1'
  []

 
  [c]
    type = ParsedMaterial
    f_name = c
    function = 'if(hSn2>=0.8, 1, 0)'
    material_property_names = 'hSn2'
  []
   
  [chSn2]
    type = ParsedMaterial
    f_name = chSn2
    function = c*hSn2
    material_property_names = 'c hSn2'
  []
 
 

  [D_original]
    type = CompositeMobilityTensor
    tensors = 'DSn1 DSn2 DCu6Sn5 DCu DGB Dint1 Dint2 DintSn DintCu'
    weights = 'ahSn1 chSn2 hGrain1 hCu hGB hint1 hint2 hintSn hintCu'
    M_name = D_original
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
   
  []

 
 
  [D]
    type = CompositeMobilityTensor
    tensors = 'D_original'
    weights = 'adjustment'
    M_name = D
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
    outputs = exodus
   
  []

  
  
  # Coefficients for diffusion equation
  [DhSn1]
    type = CompositeMobilityTensor
    tensors = 'D'
    weights = 'hSn1'
    M_name = DhSn1
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
    
  []
  [DhGrain1]
    type = CompositeMobilityTensor
    tensors = 'D'
    weights = 'hGrain1'
    M_name = DhGrain1
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
    
  []

   [DhSn2]
    type = CompositeMobilityTensor
    tensors = 'D'
    weights = 'hSn2'
    M_name = DhSn2
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
    
  []

   [DhCu]
    type = CompositeMobilityTensor
    tensors = 'D'
    weights = 'hCu'
    M_name = DhCu
    args = 'etaSn1 etaGrain1 etaSn2 etaCu'
  
  []  

  
  

  # 1/ohm-um
  [conductivity_consts]
    type = GenericConstantMaterial
    prop_names = 'condSn condCu6Sn5 condCu'
    # prop_values = '1.1e1  1.75e1     1.7e2'
   prop_values = '9.09  5.7       58.8'
  []
  [./conductivity]
    type = DerivativeParsedMaterial
    f_name = conductivity
    function = '(hGrain1*condCu6Sn5+hSn1*condSn+hSn2*condSn+hCu*condCu)'
    material_property_names = 'hSn1 hGrain1 hSn2 hCu condSn condCu6Sn5 condCu'
    outputs = exodus
    derivative_order = 2
  [../]

    [ChargeConsts]
      type = GenericConstantMaterial
      prop_names= '  kB           eCharge'
      prop_values = '1.380e-23    1.60e-19'
    []

    [EffectiveCharges]
      type = GenericConstantMaterial
      prop_names = 'ZSn      ZCu   ZCu6Sn5 ZSnEta ZCuEta'
      prop_values = '-17.97 -9.10 -13.93   -8.98  -7.87 '
      # prop_values = '-17.97 -9.10   -50   -8.98  -7.87 '
      #prop_values = '-17.97 -9.10 -20.93   -8.98  -7.87 '
      # prop_values = '17.97 9.10 13.93'
    []
  
    [Charge]
      type = ParsedMaterial
      f_name = Z
      function = 'hSn1*ZSn + hCu*ZCu + hGrain1*ZCu6Sn5 + hSn2*ZSn '#+ hSn*hGrain1*ZSnEta + hSn*hGrain2*ZSnEta 
                    #  + hGrain1*hCu*ZCuEta + hGrain2*hCu*ZCuEta'
      # function = '-50' # Tu 2005, Park 2013
      args = 'etaSn1 etaCu etaSn2 etaGrain1'
      material_property_names = 'ZSn ZCu ZCu6Sn5 ZSnEta ZCuEta hSn1(etaSn1) hGrain1(etaGrain1) hSn2(etaSn2) hCu(etaCu)'
      outputs = exodus
    []
  
    [DiffusivityWeightEM]
      type = DerivativeParsedMaterial
      f_name = D_weightEM
      # function = '1/kB/T*eCharge*(hSn*cSn*ZSn + hCu*cCu*ZCu + hGrain1*cGrain1*ZCu6Sn5 + hGrain2*cGrain2*ZCu6Sn5)'
      function = '1/kB/473*eCharge*(hSn1*(1-cSn1)*Z + hCu*(1-cCu)*Z + hGrain1*(1-cGrain1)*Z + hSn2*(1-cSn2)*Z)'
      # function = '1/kB/473*eCharge*(hSn*(0.08 + cSn)*ZSn + hCu*cCu*ZCu + hGrain1*cGrain1*ZCu6Sn5 + hGrain2*cGrain2*ZCu6Sn5)'
      args = 'cSn1 cCu cGrain1 cSn2 etaSn1 etaCu etaGrain1 etaSn2'
      material_property_names = 'kB eCharge ZSn ZCu ZCu6Sn5 hSn1(etaSn1) hCu(etaCu) hGrain1(etaGrain1) hSn2(etaSn2) Z'
      outputs = exodus
      derivative_order = 1
    []
  
    [DiffusivityEM]
      type = CompositeMobilityTensor
      tensors = 'D'
      weights = 'D_weightEM'
      args = 'cSn1 cCu cGrain1 cSn2 etaSn1 etaCu etaGrain1 etaSn2'
      M_name = D_EM
      #derivative_order = 1
    []


  [phase]
    type = ParsedMaterial
    f_name = phase
    function = 'etaCu*-2 + etaSn1 *2 + etaGrain1 * 6 + etaSn2 * 5'
    args = 'etaCu etaSn1 etaGrain1 etaSn2'
    outputs = exodus
  []

  [constraintCheck]
    type = ParsedMaterial 
    f_name = constraint
    function = 'etaCu + etaSn1 + etaGrain1 + etaSn2'
    args = 'etaCu etaSn1 etaGrain1 etaSn2'
    outputs = exodus
  []  


[]

[Kernels]

   # Poisson's Equation for Electric Field
   [poisson]
    type = MatDiffusion
    variable = psi
    diffusivity = conductivity
  []


  #Kernels for diffusion equation
  [diff_time]
    type = TimeDerivative
    variable = c
    save_in = dcdt
  []

  [diff_cSn1]
    type = MatAnisoDiffusion
    variable = c
    diffusivity = 'DhSn1'
    v = cSn1
  []

  [diff_cGrain1]
    type = MatAnisoDiffusion
    variable = c
    diffusivity = 'DhGrain1'
    v = cGrain1
  []

  [diff_cSn2]
    type = MatAnisoDiffusion
    variable = c
    diffusivity = 'DhSn2'
    v = cSn2
  []

  [diff_cCu]
    type = MatAnisoDiffusion
    variable = c
    diffusivity = 'DhCu' 
    v = cCu
  []

  [diff_EM]
    type = MatAnisoDiffusion
    variable = c
    diffusivity = D_EM
    v = psi
  []


  # Kernels for Allen-Cahn equation for etaSn

# Kernels for Allen-Cahn equation for etaSn

  [detaS1dt]
    type = TimeDerivative
    variable = etaSn1
    save_in = detaSn1dt
  []
  
  [ACBulkFSn] # h'f + w*g
    type = KKSMultiACBulkF
    variable = etaSn1
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gSn1
    eta_i = etaSn1
    wi = 4e2
    args = 'cSn1 cGrain1 cSn2 cCu etaGrain1 etaSn2 etaCu psi'
    save_in = acbulkf
  []
  [ACBulkCSn] # h*dfdc * c
    type = KKSMultiACBulkC
    variable = etaSn1
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaSn1
    args = 'etaGrain1 etaSn2 etaCu psi'
    save_in = acbulkc
  []

  #[ACInterfaceSn1] # eps_sq *laplacian(phi)
    #type = ACInterface
    #variable = etaSn1
    #kappa_name = kappa
    #save_in = acinterface
  #[]

  [ACInterfaceSn1_Grain1] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn1
    kappa_name = kappa
    coupled_variables = etaGrain1
    save_in = acinterface
  []

  [ACInterfaceSn1_Sn2] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn1
    kappa_name = kappa
    coupled_variables = etaSn2
    save_in = acinterface
  []
  [ACInterfaceSn1_Cu] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn1
    kappa_name = kappa
    coupled_variables = etaCu
    save_in = acinterface
  []
  
  
  [multiplerSn] # lambda from variational derivative
    type = MatReaction
    variable = etaSn1
    v = lambda
    mob_name = L
    save_in = lambdakernel
  []

  # Kernels for Allen-Cahn equation for Grain 1
  [detaGrain1dt]
    type = TimeDerivative
    variable = etaGrain1
  []
  [ACBulkFGrain1]
    type = KKSMultiACBulkF
    variable = etaGrain1
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gGrain1
    eta_i = etaGrain1
    wi = 4e2
    args = 'cSn1 cGrain1 cSn2 cCu etaSn1 etaSn2 etaCu psi'
  []
  [ACBulkCGrain1]
    type = KKSMultiACBulkC
    variable = etaGrain1
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaGrain1
    args = 'etaSn1 etaSn2 etaCu psi'
  []

  
  [ACInterfaceGrain1_Sn1]
    type = ACInterface
    variable = etaGrain1
    kappa_name = kappa
    coupled_variables = etaSn1
  []
  [ACInterfaceGrain1_Sn2]
    type = ACInterface
    variable = etaGrain1
    kappa_name = kappa
    coupled_variables = etaSn2
  []

   [ACInterfaceGrain1_Cu]
    type = ACInterface
    variable = etaGrain1
    kappa_name = kappa
    coupled_variables = etaCu
  []
 
  
  [multiplerGrain1]
    type = MatReaction
    variable = etaGrain1
    v = lambda
    mob_name = L
  []

  # Kernels for Allen-Cahn equation for Sn 2
  [detaSn2dt]
    type = TimeDerivative
    variable = etaSn2
  []
  [ACBulkFSn2]
    type = KKSMultiACBulkF
    variable = etaSn2
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gSn2
    eta_i = etaSn2
    wi = 4e2
    args = 'cSn1 cGrain1 cSn2 cCu etaSn1 etaGrain1 etaCu psi'
  []

  [ACBulkCSn2]
    type = KKSMultiACBulkC
    variable = etaSn2
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaSn2
    args = 'etaSn1 etaGrain1 etaCu psi'
  []

  
 
  [ACInterfaceSn2_Grain1] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn2
    kappa_name = kappa
    save_in = acinterface
    coupled_variables = etaGrain1
  []
  [ACInterfaceSn2_Sn1] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn2
    kappa_name = kappa
    coupled_variables = etaSn1
  []

  [ACInterfaceSn2_Cu] # eps_sq *laplacian(phi)
    type = ACInterface
    variable = etaSn2
    kappa_name = kappa
    coupled_variables = etaCu
  []


  [multiplerSn2]
    type = MatReaction
    variable = etaSn2
    v = lambda
    mob_name = L
  []


  # Kernels for the Lagrange multiplier equation
  [mult_lambda]
    type = MatReaction
    variable = lambda
    mob_name = 4
  []
  [mult_ACBulkF_Sn1]
    type = KKSMultiACBulkF
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gSn1
    eta_i = etaSn1
    wi = 4e2
    mob_name = 1
    args = 'cSn1 cGrain1 cSn2 cCu etaGrain1 etaSn2 etaCu psi'
  []
  [mult_ACBulkC_Sn1]
    type = KKSMultiACBulkC
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaSn1
    args = 'etaGrain1 etaSn2 etaCu psi'
    mob_name = 1
  []
  [mult_CoupledACint_Sn1]
    type = SimpleCoupledACInterface
    variable = lambda
    v = etaSn1
    kappa_name = kappa
    mob_name = 1
  []

  [mult_ACBulkF_Grain1]
    type = KKSMultiACBulkF
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gGrain1
    eta_i = etaGrain1
    wi = 4e2
    mob_name = 1
    args = 'cSn1 cGrain1 cSn2 cCu etaSn1 etaSn2 etaCu psi'
  []
  [mult_ACBulkC_Grain1]
    type = KKSMultiACBulkC
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaGrain1
    args = 'etaSn1 etaSn2 etaCu psi'
    mob_name = 1
  []
  [mult_CoupledACint_Grain1]
    type = SimpleCoupledACInterface
    variable = lambda
    v = etaGrain1
    kappa_name = kappa
    mob_name = 1
  []

  ##
    [mult_ACBulkF_Sn2]
    type = KKSMultiACBulkF
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gSn2
    eta_i = etaSn2
    wi = 4e2
    mob_name = 1
    args = 'cSn1 cGrain1 cSn2 cCu etaSn1 etaGrain1 etaCu psi'
  []
  [mult_ACBulkC_Sn2]
    type = KKSMultiACBulkC
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaSn2
    args = 'etaSn1 etaGrain1 etaCu psi'
    mob_name = 1
  []
  [mult_CoupledACint_Sn2]
    type = SimpleCoupledACInterface
    variable = lambda
    v = etaSn2
    kappa_name = kappa
    mob_name = 1
  []

 ##
   [mult_ACBulkF_Cu]
    type = KKSMultiACBulkF
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    gi_name = gCu
    eta_i = etaCu
    wi = 4e2
    mob_name = 1
    args = 'cSn1 cGrain1 cSn2 cCu etaSn1 etaGrain1 etaSn2 psi'
  []
  [mult_ACBulkC_Cu]
    type = KKSMultiACBulkC
    variable = lambda
    Fj_names = 'FSn1 FGrain1 FSn2 FCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    cj_names = 'cSn1 cGrain1 cSn2 cCu'
    eta_i = etaCu
    args = 'etaSn1 etaGrain1 etaSn2 psi'
    mob_name = 1
  []
  [mult_CoupledACint_Cu]
    type = SimpleCoupledACInterface
    variable = lambda
    v = etaCu
    kappa_name = kappa
    mob_name = 1
  []


  # Kernels for constraint equation eta1 + eta2 + eta3 = 1
  # eta3 is the nonlinear variable for the constraint equation
  [etaCuReaction]
    type = MatReaction
    variable = etaCu
    mob_name = 1
  []
  [etaSn1reaction]
    type = MatReaction
    variable = etaCu
    v = etaSn1
    mob_name = 1
  []
  [etaGrain1reaction]
    type = MatReaction
    variable = etaCu
    v = etaGrain1
    mob_name = 1
  []
  [etaSn2reaction]
    type = MatReaction
    variable = etaCu
    v = etaSn2
    mob_name = 1
  []
  [one]
    type = BodyForce
    variable = etaCu
    value = -1.0
  []

  # Phase concentration constraints
  [chempot12]
    type = KKSPhaseChemicalPotential
    variable = cSn1
    cb = cGrain1
    fa_name = FSn1
    fb_name = FGrain1
  []
  [chempot23]
    type = KKSPhaseChemicalPotential
    variable = cGrain1
    cb = cSn2
    fa_name = FGrain1
    fb_name = FSn2
  []
  [chempot34]
    type = KKSPhaseChemicalPotential
    variable = cSn2
    cb = cCu
    fa_name = FSn2
    fb_name = FCu
  []


  [phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = cCu
    cj = 'cSn1 cGrain1 cSn2 cCu'
    hj_names = 'hSn1 hGrain1 hSn2 hCu'
    etas = 'etaSn1 etaGrain1 etaSn2 etaCu'
    c = c
  []

[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 40
  nl_max_its = 30
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-5
  nl_abs_tol = 1.0e-7

  # automatic_scaling = true
  # compute_scaling_once = false

  # Start off with this small timestep to resolve sharp interfaces.
  # num_steps = 20
  # dt = 1e-5

  num_steps = 120000
  dt = 1e-6

[]

[Preconditioning]
  active = 'full'
  #active = 'mydebug'
  #active = ''
  [full]
    type = SMP
    full = true
  []
  # [mydebug]
  #   type = FDP
  #   full = true
  # []
[]

[Outputs]
  exodus = true
  checkpoint = true
  interval = 100
[]

[Debug]
  show_var_residual_norms = true
[]



