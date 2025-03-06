model SolidOxideFuelCell
  extends Modelica.Icons.Package;

  model Examples
    extends Modelica.Icons.ExamplesPackage;

    model SOFC_Validation
      extends Modelica.Icons.Example;
  Components.SOFC FuelCell annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
    equation

    end SOFC_Validation;
  equation

  end Examples;

  model Components
    extends Modelica.Icons.Package;

    model SOFC
    
      replaceable parameter 
        SolidOxideFuelCell.Records.Stack_5kW SOFC_Data constrainedby 
        SolidOxideFuelCell.Records.PartialRecord
        "Specific SOFC data record" annotation (
        choicesAllMatching=true,
        Dialog(group="Fuel cell model selection"),
        Placement(transformation(extent={{-10,-10},{10,10}})));
    
    
      import Modelica.Units.SI;
    //Constants
      constant SI.FaradayConstant F = Modelica.Constants.F                          "Faraday constant [C/mol]";
      constant Real R = Modelica.Constants.R                                        "Molar gas constant [J/(mol.K)]";
      constant Real sigma = Modelica.Constants.sigma                                "Stefan-Boltzmann constant";
      constant SI.MolarMass MM_air = 28.95*1e-3                                     "Molar mass air [kg/mol]";
      constant SI.MolarMass MM_H2 = 2e-3                                            "Molar mass hydrogen  [kg/mol]";
      constant SI.MolarMass MM_H2O = 18e-3                                          "Molar mass water [kg/mol]";
      constant SI.MolarMass MM_O2 = 31.999*1e-3                                     "Molar Mass of oxygen [kg/mol]";
      constant SI.MolarMass MM_N2 = 28.01340*1e-3                                   "Molar Mass of nitrogen [kg/mol]";
      constant SI.MolarEnthalpy h_molar_H2 = 285.83e3                               "Molar specific heating value hydrogen [J/mol]";
      constant SI.SpecificEnergy h_inf_H2= 1.19988e+8                               "Gross calorific heating value hydrogen [J/kg]";
    //Input variable
      SI.Current I_fc                                                               "SOFC current [A]";
    //Combi time table for dynamic input current signal
      import Modelica.Blocks.Sources.CombiTimeTable;
      Modelica.Blocks.Sources.CombiTimeTable I_signal(
       table = [ 0, 0;
                600, 0;
                600, 100;
                7200, 100;
                7200, 30;
                8000, 30],
      columns = {1,2}, extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints);
    //FC specifications
      parameter Integer N_cells = SOFC_Data.N_cells                                 "Number of cells" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Volume V_a = SOFC_Data.V_a                                       "Volume anode gas channel [m^3]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Volume V_c = SOFC_Data.V_c                                       "Volume anode gas channel [m^3]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Length delta = SOFC_Data.delta                                   "Length cell [m]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Area  A_electrolyte = SOFC_Data.A_electrolyte                    "Electrolyte area [m^2]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Area A_interconnect = SOFC_Data.A_interconnect                   "Interconnect area [m^2]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Current I_limit = SOFC_Data.I_limit                              "Limiting Current [A]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Temperature T_op = SOFC_Data.T_op                                "FC operating temperautre [K]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications", defautl = SOFC_Data.T_op));
      parameter SI.Area A_ast_out = SOFC_Data.A_ast_out                             "Outer area air supply tube [m^2]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Area A_ast_in = SOFC_Data.A_ast_in                               "Inner area air supply tube [m^2]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Volume V_ast_cave = SOFC_Data.V_ast_cave                         "Volume air supply channel cave [m^3]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Mass M_ast = SOFC_Data.M_ast                                     "Mass air supply tube [kg]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Area A_cell = SOFC_Data.A_cell                                   "Area cell [m^2]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Volume V_ann = SOFC_Data.V_ann                                   "Volume of annular space between ast and cell [m^3]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.Mass M_cell = SOFC_Data.M_cell                                   "Mass single cell[kg]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
      parameter SI.SpecificHeatCapacity c_p_cell = SOFC_Data.c_p_cell               "Specific heat capacity single cell[J/(kgK)]" annotation(dialog(group = "Data sheet", tab = "Fuel cell specifications"));
    //Tunable parameters
      //parameter Real C_eq = 0.0004;
      parameter SI.Voltage xi_0 = SOFC_Data.xi_0                                    "Tunable parameter [V]" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real xi_1 = SOFC_Data.xi_1                                          "Tunable parameter [V/K]" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real k_e = SOFC_Data.k_e                                            "Tunable parameter [V/K]" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real a_electrolyte = SOFC_Data.a_electrolyte                        "Tunable parameter" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real b_electrolyte = SOFC_Data.b_electrolyte                        "Tunable parameter" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real a_interconnect = SOFC_Data.a_interconnect                      "Tunable parameter" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real b_interconnect = SOFC_Data.b_interconnect                      "Tunable parameter" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real k_a1 = SOFC_Data.k_a1                                          "Tunable parameter [A/K]" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
      parameter Real k_a2 = SOFC_Data.k_a2                                          "Tunable parameter" annotation(dialog(group = "Tunable parameters", tab = "Fuel cell specifications"));
    //Operaring parameters
      parameter SI.Pressure p_a = SOFC_Data.p_a                                     "Anode pressure [Pa]" annotation(dialog(group = "Operating parameters"));
      parameter SI.Pressure p_c = SOFC_Data.p_c                                     "Cathode pressure [Pa]"annotation(dialog(group = "Operating parameters"));
      parameter SI.Temperature T_air_in = SOFC_Data.T_air_in                        "Temperature air inlet [K]" annotation(dialog(group = "Operating parameters"));
      parameter SI.Temperature T_fuel_in = SOFC_Data.T_fuel_in                      "Temperature fuel inlet [K]" annotation(dialog(group = "Operating parameters"));
      parameter SI.Temperature T_amb = SOFC_Data.T_amb                              "Ambient temperature [K]" annotation(dialog(group = "Operating parameters"));
      parameter SI.MolarFlowRate n_H2_in_set = SOFC_Data.n_H2_in                    "Inlet molar flow hydrogen [mol/s]" annotation(dialog(group = "Constand fuel flow operation"));
      parameter SI.MolarFlowRate n_H2O_in_set = SOFC_Data.n_H2O_in                  "Inlet molar flow water [mol/s]" annotation(dialog(group = "Constand fuel flow operation"));
      parameter SI.MolarFlowRate n_a_set = SOFC_Data.n_a                            "Anode molar flow (fuel flow) [mol/s]" annotation(dialog(group = "Constand fuel flow operation"));
      parameter SI.MolarFlowRate n_c_set = SOFC_Data.n_c                            "Cathode flow (Air flow) [mol/s]" annotation(dialog(group = "Constand fuel flow operation"));
      //Constant fuel utilization operation
      parameter Real u_fuel_set = SOFC_Data.u_fuel                                  "Fuel utilization rate as a decimal number" annotation(dialog(group = "Constand fuel utiilization operation"));
    //Pneumatical submodel
      //Electrodes
      SI.MolarFlowRate n_a                                                          "Molar flow rate anode inlet [mol/s]";
      SI.MolarFlowRate n_c                                                          "Molar flow rate cathode inlet [mol/s]";
    //Fuel utilization
      Real u_fuel                                                                   "Amount of hydrogen in inlet flow rate which takes part in reaction [%]";
    //Hydrogen
      SI.MolarFlowRate n_H2_in                                                      "Molar flow rate hydrogen inlet [mol/s]";
      SI.MassFlowRate m_H2_in                                                       "Mass flow rate hydrogen inlet [kg/s]";
      Real x_H2_in                                                                  "Molar fraction H2 inlet";
      SI.Pressure p_H2_in                                                           "Inlet partial pressure hydrogen [Pa]";
      SI.Pressure p_H2_ch(start=p_H2_start)                                         "Partial pressure hydrogen in gas channel";
      SI.Pressure p_H2_eff                                                          "Effective partial pressure hydrogen [Pa]";
      SI.Pressure p_H2_out                                                          "Effective partial pressure hydrogen outlet [Pa]";
      Real x_H2_out                                                                 "Molar fraction hydrogen outlet";
      SI.MolarFlowRate n_H2_out                                                     "Molar flow hydrogen outlet [mol/s]";
      SI.MassFlowRate m_H2_out                                                      "Mass flow rate hydrogen outlet [kg/s]";
      SI.MolarFlowRate n_H2_reac                                                    "Molar flow hydrogen involved in reaction [mol/s]";
    //Water
      SI.MolarFlowRate n_H2O_in                                                     "Molar flow rate water inlet [mol/s]";
      SI.MassFlowRate m_H2O_in                                                      "Mass flow rate water inlet [kg/s]";
      Real x_H2O_in                                                                 "Molar fraction H2O inlet";
      SI.Pressure p_H2O_in                                                          "Inlet partial pressure water [Pa]";
      SI.Pressure p_H2O_ch(start=p_H2O_start)                                       "Partial pressure water in gas channel";
      SI.Pressure p_H2O_eff                                                         "Effective partial pressure water [Pa]";
      SI.Pressure p_H2O_out                                                         "Effective partial pressure water outlet [Pa]";
      Real x_H2O_out                                                                "Molar fraction water outlet";
      SI.MolarFlowRate n_H2O_out                                                    "Molar flow water outlet";
      SI.MassFlowRate m_H2O_out                                                     "Mass flow rate water outlet [kg/s]";
      SI.MolarFlowRate n_H2O_reac                                                   "Molar flow water due to reaction [mol/s]";
    //Oxygen
      SI.MassFlowRate m_O2_in                                                       "Mass flow rate oxygen inlet [kg/s]";
      Real x_O2_in                                                                  "Molar fraction oxygen inlet";
      SI.Pressure p_O2_in                                                           "Inlet partial pressure oxygen [Pa]";
      SI.MolarFlowRate n_O2_in                                                      "Molar flow oxygen inlet";
      SI.Pressure p_O2_ch(start=p_O2_start)                                         "Partial pressure oxygen in gas channel";
      SI.Pressure p_O2_eff                                                          "Effective partial pressure oxygen [Pa]";
      SI.Pressure p_O2_out                                                          "Effective partial pressure oxygen outlet [Pa]";
      Real x_O2_out                                                                 "Molar fraction oxygen outlet";
      SI.MolarFlowRate n_O2_out                                                     "Molar flow oxygen outlet";
      SI.MassFlowRate m_O2_out                                                      "Mass flow rate oxygen outlet [kg/s]";
      SI.MolarFlowRate n_O2_reac                                                    "Molar flow oxygen involved in reaction [mol/s]";
    //Nitrogen
      SI.MassFlowRate m_N2_in                                                       "Mass flow rate nitrogen inlet [kg/s]";
      Real x_N2_in                                                                  "Molar fraction nitrogen inlet";
      SI.Pressure p_N2_in                                                           "Inlet partial pressure nitrogen [Pa]";
      SI.MolarFlowRate n_N2_in                                                      "Molar flow nitrogen inlet";
      SI.Pressure p_N2_out                                                          "Effective partial pressure oxygen outlet [Pa]";
      Real x_N2_out                                                                 "Molar fraction oxygen outlet";
      SI.MolarFlowRate n_N2_out                                                     "Molar flow oxygen outlet";
      SI.MassFlowRate m_N2_out                                                      "Mass flow rate oxygen outlet [kg/s]";
    //Electrochemical submodel
      //Ohmic loss
      SI.Resistance R_electrolyte                                                   "Ohmic resistance electrolyte [Ohm]";
      SI.Voltage U_ohm_electrolyte                                                  "Ohmic voltage drop electrolyte [Ohm]";
      SI.Resistance R_interconnect                                                  "Ohmic resistance interconnect [Ohm]";
      SI.Voltage U_ohm_interconnect                                                 "Ohmic voltage drop interconnect [Ohm]";
      SI.Voltage U_ohm                                                              "Ohmic voltage drop cell [Ohm]";
    //Concentration loss
      SI.Voltage U_conc                                                             "Concentration voltage drop cell [V]";
      SI.Resistance R_conc                                                          "Concentration resistance cell [Ohm]";
    //Activation loss
      SI.Current i_0                                                                "Exchange current [A]";
      SI.Voltage U_act1                                                             "Temperature and current dependant activation voltage drop [V]";
      SI.Voltage U_act0                                                             "Temperature and not current dependandt activation voltage drop [V]";
      SI.Voltage U_act                                                              "Activation voltage drop cell [V]";
      SI.Resistance R_act                                                           "Activation resistance cell [Ohm]";
    //Reversible voltage
      SI.Voltage E_cell                                                             "Reversible voltage cell [V]";
      SI.Voltage E_0                                                                "Reversible voltage cell at standard state [V]";
    //Output voltage
      SI.Voltage U_cell                                                             "Cell voltage [V]";
    //Double layer charging effect
      SI.Voltage U_dl                                                               "Voltage drop due to double layer charging effect [V]";
    //Thermal submodel
      //Air supply tube
      //Radiation air supply tube
      SI.Temperature T_ast(start=T_air_in, displayUnit="K")                         "Air supply tube temperature [K]";
      SI.HeatFlowRate Q_rad_ast                                                     "Radiation heat flow [W]";
      Real eps_m_ast                                                                "Emissivity air supply tube";
      Real eps_m_ast_eff                                                            "Effective air supply tube";
    //Convection outer sourface air supply tube
      SI.HeatFlowRate Q_out_conv_ast                                                "Convection heat loss outside [W]";
      parameter Real NusseltNumber_ast_out = 4.076                                  "Nusselt number for air supply tube on the outside" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      parameter SI.Length D_h_ast_out = 0.626e-2                                    "Hydrauclic dameter air supply tube outsie [m]" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      Real alpha_ast_out                                                            "Alpha/Heat transfer coefficient cell fuel [W/(m^2K)]";
      Real k_f_air_out                                                              "Thermal conductivity air out [W/(mK)]";
    //Convection inner surface air supply tube
      SI.HeatFlowRate Q_in_conv_ast                                                 "Convection heat loss inside [W]";
      SI.Temperature T_air_ast(start=T_air_in, displayUnit="K")                     "Air temperature cell [K]";
      parameter Real NusseltNumber_ast_in = 4.364                                   "Nusselt number for air supply tube on the inside" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      parameter SI.Length D_h_ast_in = 0.396e-2                                     "Hydrauclic dameter air supply tube intsie [m]" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      Real alpha_ast_in                                                             "Alpha/Heat transfer coefficient inner surface air supply tube [W/(m^2K)]";
      Real k_f_air_in                                                               "Thermal conductivity air in [W/(mK)]";
    //Heat rate air flow in air supply tube
      SI.HeatFlowRate Q_air_ast                                                     "Heat flow air [W]";
      SI.MassFlowRate m_air                                                         "Molar flow rate air [kg/s]";
      SI.SpecificHeatCapacity c_p_air_ast                                           "Specific heat capacity ais [J/(kgK)]";
    //Energy balacne air supply tube
      SI.HeatFlowRate Q_ast                                                         "Heat rate supply tube [W]";
    //Temperature in air supply tube
      SI.SpecificHeatCapacity c_p_ast                                               "Specific heat capacity air supply tube [J/(kgK)]";
    //Cell
      //Heat rate fuel flow
      SI.HeatFlowRate Q_fuel                                                        "Heat flow fuel [W]";
      SI.ThermalConductance C_th_fuel_tot                                           "Thermal mass flow fuel in anode gas channel [W/K]";
      SI.SpecificHeatCapacity c_p_H2                                                "Specific heat capacity hydrogen cell [J/(kgK)]";
      SI.MassFlowRate m_H2_tot                                                      "Molar flow rate hydrogen cell [kg/s]";
      SI.MolarFlowRate n_H2_tot                                                     "Molar flow rate hydrogen cell [mol/s]";
      SI.ThermalConductance C_th_H2_tot                                             "Thermal mass flow rate hydrogen cell [W/kg]";
      SI.SpecificHeatCapacity c_p_H2O                                               "Specific heat capacity water cell [J/(kgK)]";
      SI.MolarFlowRate n_H2O_tot                                                    "Molar flow rate water cell [mol/s]";
      SI.MassFlowRate m_H2O_tot                                                     "Mass flow rate water cell [kg/s]";
      SI.ThermalConductance C_th_H2O_tot                                            "Thermal mass flow rate water cell [W/kg]";
    //Convection cell to fuel
      SI.HeatFlowRate Q_conv_cell_fuel                                              "Convection heat rate cell fuel [W]";
      Real alpha_cell_fuel                                                          "Alpha/Heat transfer coefficient cell fuel [W/(m^2K)]";
      Real lambda_fuel                                                              "Thermal conductivity fuel [W/(mK)]";
      Real K_f_H2                                                                   "Coefficient for thermal conductivity hydrogen cell [W/(mK)]";
      Real K_f_H2O                                                                  "coefficient for thermal conductivity water cell [W/(mK)]";
      Real lambda_H2                                                                "Thermal conductivity hydrogen cell [W/(mK)]";
      Real lambda_H2O                                                               "Thermal conductivity water cell [W/(mK)]";
      parameter Real NusseltNumber_cell_fuel = 5.9                                  "Nusselt number cell fuel" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      parameter SI.Length D_h_cell_fuel = 0.644e-2                                  "Hydrauclic dameter cell fuel side/anode side [m]" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
    //Heat rate air flow in cell
      SI.HeatFlowRate Q_air_cell                                                    "Heat rate air flow in cell [W]";
      SI.SpecificHeatCapacity c_p_air_cell                                          "Specific heat capacity air cell [J/(kgK)]";
    //Convection cell to air
      SI.HeatFlowRate Q_conv_cell_air                                               "Convection heat rate cell air [W]";
      Real alpha_cell_air                                                           "Alpha/Heat transfer coefficient cell air [W/(m^2K)]";
      Real lambda_air_cell                                                          "Thermal Conductivity air cell [W/(mK)]";
      parameter Real NusseltNumber_cell_air = 4.169/(1+0.2042)                      "Nusselt number cell air" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
      parameter SI.Length D_h_cell_air = 0.626e-2                                   "Hydrauclic dameter cell air side/cathode [m]" annotation(dialog(group = "Thermal parameters", tab = "Fuel cell specifications"));
    //Temperature air in cell
      SI.Temperature T_air_cell(start=T_air_in, displayUnit="K")                    "Temperature air cell [K]";
      SI.SpecificHeatCapacity c_p_air_ann                                           "Specific heat capacity air annual space [J/(kgK)]";
      SI.AmountOfSubstance N_air_ann                                                "Amount of substance air in cathode gas channel [mol]";
      SI.Mass M_air_ann                                                             "Mass air in annual space [kg]";
      SI.HeatCapacity C_th_air_ann                                                  "Thermal mass air in annual space [J/kg]";
    //Temperature fuel in cell
      SI.Temperature T_fuel_out(start=T_fuel_in, displayUnit="K")                       "Temperature fuel [K]";
    //Theroetical maximum chemical reaction heat
      SI.HeatFlowRate Q_max_reac                                                    "Theoretical maxium reaction heat supplied by hydrogen inlet flow [W]";
      SI.MolarEnthalpy deltah_H2                                                    "Gross caloric specific heating valur hydrogen [J/mol]";
      parameter SI.MolarEnthalpy deltah_H2_0 = 241820                               "Molar lower heating value hydrogen [J/(molK)]";
      SI.MolarEnthalpy delta_deltah_H2                                              "Temperature dependant change of gross caloric specific heating valur hydrogen [J/mol]";
      parameter Real A_H2 = 33.066178                                               "Shomate coefficient hydrogen";
      parameter Real B_H2 = -11.363417                                              "Shomate coefficient hydrogen";
      parameter Real C_H2 = 11.432816                                               "Shomate coefficient hydrogen";
      parameter Real D_H2 = -2.772874                                               "Shomate coefficient hydrogen";
      parameter Real E_H2 = -0.158558                                               "Shomate coefficient hydrogen";
      parameter Real A_H2O = 30.09200                                               "Shomate coefficient water";
      parameter Real B_H2O = 6.832514                                               "Shomate coefficient water";
      parameter Real C_H2O = 6.793435                                               "Shomate coefficient water";
      parameter Real D_H2O = -2.534480                                              "Shomate coefficient water";
      parameter Real E_H2O = 0.082139                                               "Shomate coefficient water";
      parameter SI.Temperature T_0 = 25+273.15                                      "Standard temperature [K]";
      parameter Integer N_int = 10                                                  "Number of intervals";
      Real deltaT                                                                   "Stepsize for iteration";
      Real T_start[N_int]                                                           "Startdtemperatur for specific integration [K]";
      Real T_end[N_int]                                                             "Endtemperatur for specific integration [K]";
      Real Delta_c_p_H2_start[N_int]                                                "Specific heat capacitie at end of specific integral [kJ/(kgK)]";
      Real Delta_c_p_H2_end[N_int]                                                  "Specific heat capacitie at end of specific integral [kJ/(kgK)]";
    //Electrical power fuel cell
      SI.Power P_cell                                                               "Electrical power output fuel cell single cell [W]";
    //Thermal mass cell
      SI.HeatCapacity C_th_cell                                                     "Thermal mass [J/kg]";
    //Energy balance
      SI.HeatFlowRate Q_cell                                                        "Heat rate cell [W]";
    //Fuel cell temperatures
      SI.Temperature T_cell(start=T_fuel_in, displayUnit="K")                       "Cell temperature [K]";
      SI.HeatCapacity C_th_fuel                                                     "Thermal mass fuel in anode gas channel [J/K]";
      SI.HeatCapacity C_th_H2                                                       "Thermal mass hydrogen in anode gas channel [J/K]";
      SI.HeatCapacity C_th_H2O                                                      "Thermal mass water in anode gas channel[W/K]";
      SI.AmountOfSubstance N_H2                                                     "Amount of substance hydrogen in anode gas channel [mol]";
      SI.Mass M_H2                                                                  "Mass hydrogen in anode gas channel [kg]";
      SI.AmountOfSubstance N_H2O                                                    "Amount of substance water in anode gas channel [mol]";
      SI.Mass M_H2O                                                                 "Mass water in anode gas channel [kg]";
    //Stack variables
      SI.Temperature T_stack(displayUnit="K")                                       "Fuel cell stack temperature [K]";
      SI.Power P_stack                                                              "Electrical power output fuel cell stack [W]";
      SI.Voltage U_stack                                                            "Stack voltage [V]";
      SI.Voltage E_stack                                                            "Reversible voltage stack [V]";
      SI.MolarFlowRate n_H2_in_stack                                                "Molar inlet flow stack hydrogen stack [mol/s]";
      SI.MolarFlowRate n_H2O_in_stack                                               "Molar inlet flow stack water stack [mol/s]";
      SI.MolarFlowRate n_O2_in_stack                                                "Molar inlet flow stack oxygen stack[mol/s]";
      SI.MassFlowRate m_H2_in_stack                                                 "Mass inlet flow stack hydrogen stack [kg/s]";
      SI.MassFlowRate m_H2O_in_stack                                                "Mass inlet flow stack water stack [kg/s]";
      SI.MassFlowRate m_O2_in_stack                                                 "Mass inlet flow stack oxygen stack [kg/s]";
      SI.MolarFlowRate n_H2_out_stack                                               "Molar outlet flow stack hydrogen stack [mol/s]";
      SI.MolarFlowRate n_H2O_out_stack                                              "Molar outlet flow stack water stack [mol/s]";
      SI.MolarFlowRate n_O2_out_stack                                               "Molar outlet flow stack oxygen stack [mol/s]";
      SI.MassFlowRate m_H2_out_stack                                                "Mass outlet flow stack hydrogen stack [kg/s]";
      SI.MassFlowRate m_H2O_out_stack                                               "Mass outlet flow stack water stack [kg/s]";
      SI.MassFlowRate m_O2_out_stack                                                "Mass outlet flow stack oxygen stack [kg/s]";
      Real eta_far                                                                  "Faraday efficiency [%] ";
      Real eta_el                                                                   "Electrical efficiency [%]";
    //Start variables for integration
      parameter SI.Current I_init = 0                                               "Auxiliary variable [A]";
      parameter SI.Pressure p_H2_start = 101325 ;      //p_H2_in - I_init*p_a/(4*n_a*F)            "Initial value for p_H2_ch [Pa]";
      parameter SI.Pressure p_H2O_start = 101325 ;//p_H2O_in + I_init*p_a/(4*n_a*F)          "Initial value for p_H2O_ch [Pa]";
      parameter SI.Pressure p_O2_start = 101325;    
//p_O2_in - I_init*p_c/(8*n_c*F)             "Initial value for p_H2O_ch [Pa]";
      Modelica.Blocks.Sources.Ramp ramp(
        height=140,
        duration=100,
        offset=0.1,
        startTime=0)
        annotation (Placement(transformation(extent={{-10,32},{10,52}})));
    
    algorithm
    
     //Calculation of the temperature dependant heating value of hydrogen
//Stepsize
      deltaT := (T_cell - T_0)/N_int;
//Initialization
      delta_deltah_H2 := 0.0;
    
      for i in 1:N_int loop
//Start and end conditions for each interval
        T_start[i] := T_0 + (i - 1)*deltaT;
        T_end[i] := T_0 + i*deltaT;
//Specific heat capacities with Shomate equations
        Delta_c_p_H2_start[i] := (A_H2O + B_H2O * (T_start[i]/1000) +C_H2O*(T_start[i]/1000)^2 +D_H2O*(T_start[i]/1000)^3 +E_H2O/(T_start[i]/1000)^2) - (A_H2 + B_H2*(T_start[i]/1000) + C_H2*(T_start[i]/1000)^2 + D_H2*(T_start[i]/1000)^3 + E_H2/(T_start[i]/1000)^2);
        Delta_c_p_H2_end[i] := (A_H2O + B_H2O * (T_end[i]/1000) + C_H2O*(T_end[i]/1000)^2 +D_H2O*(T_end[i]/1000)^3 +E_H2O/(T_end[i]/1000)^2) -(A_H2 + B_H2*(T_end[i]/1000) +C_H2*(T_end[i]/1000)^2 +D_H2*(T_end[i]/1000)^3 + E_H2/(T_end[i]/1000)^2);
//Heating value with specific integration
        delta_deltah_H2 := delta_deltah_H2 + deltaT/2*(Delta_c_p_H2_start[i] + Delta_c_p_H2_end[i]);
      end for;
    equation
    
        I_fc = I_signal.y[2];
        n_H2_in = n_H2_in_set;
        n_H2O_in = n_H2O_in_set;
        n_a = n_a_set;
        n_c = n_c_set;
        u_fuel = n_H2_reac/n_H2_in;
//Pneumatical submodel
//Hydrogen
      m_H2_in = n_H2_in*MM_H2;
      if n_a <= 0 then
        x_H2_in = 0;
      else
        x_H2_in = n_H2_in/n_a;
      end if;
      p_H2_in = p_a*x_H2_in;
      der(p_H2_ch) = 2*(p_H2_in - p_H2_ch)*(n_a*R*T_cell/(V_a*p_a)) - I_fc*R*T_cell/(V_a*2*F);
      p_H2_eff = p_H2_ch*(1-I_fc/I_limit);
      p_H2_out = 2*p_H2_ch - p_H2_in;
      n_H2_out = n_H2_in - n_H2_reac;
      n_H2_reac = I_fc/(2*F);
      m_H2_out = n_H2_out*MM_H2;
      x_H2_out = n_H2_out/(n_H2_out + n_H2O_out);
//Water
      m_H2O_in = n_H2O_in*MM_H2O;
      if n_a <= 0 then
        x_H2O_in = 0;
      else
        x_H2O_in = n_H2O_in/n_a;
      end if;
      p_H2O_in = p_a*x_H2O_in;
      der(p_H2O_ch) = 2*(p_H2O_in - p_H2O_ch)*(n_a*R*T_cell/(V_a*p_a)) + I_fc*R*T_cell/(V_a*2*F);
      p_H2O_eff = p_H2O_ch + p_H2_ch*(I_fc/I_limit);
      p_H2O_out = 2*p_H2O_ch - p_H2O_in;
      n_H2O_out = n_H2O_in + n_H2O_reac;
      m_H2O_out = n_H2O_out*MM_H2O;
      n_H2O_reac = I_fc/(2*F);
      x_H2O_out = n_H2O_out/(n_H2O_out + n_H2O_out);
//Oxygen
      m_O2_in = n_O2_in*MM_O2;
      x_O2_in = 0.22;
      p_O2_in = x_O2_in*p_c;
      n_O2_in = x_O2_in*n_c;
      der(p_O2_ch) = 2*(p_O2_in - p_O2_ch)*(n_c*R*T_cell/(V_c*p_c)) - I_fc*R*T_cell/(V_c*4*F);
      p_O2_eff = (exp(I_fc/I_limit*1/1.1*0.2485)*p_O2_ch) - (((exp(I_fc/I_limit*1/1.1*0.2485)) - 1)*p_c);
      p_O2_out = 2*p_O2_ch - p_O2_in;
      n_O2_out = n_O2_in - n_O2_reac;
      m_O2_out = n_O2_out*MM_O2;
      n_O2_reac = I_fc/(4*F);
      x_O2_out = n_O2_out/(n_N2_out + n_O2_out);
//Nitrogen
      m_N2_in = n_N2_in*MM_N2;
      x_N2_in = 1 - x_O2_in;
      p_N2_in = x_N2_in*p_c;
      n_N2_in = x_N2_in*n_c;
      n_N2_out = n_N2_in;
      x_N2_out = n_N2_out/(n_N2_out + n_O2_out);
      p_N2_out = p_c*x_N2_out;
      m_N2_out = n_N2_out*MM_N2;
//Electrochemical submodel
//Ohmic loss
      R_electrolyte = 1.39*10*a_electrolyte*exp(10350/T_cell/1.5)*delta/A_electrolyte;
      U_ohm_electrolyte = R_electrolyte*I_fc;
      R_interconnect = (1.39*10*a_interconnect*exp(b_interconnect/T_cell/1.5))*delta/A_interconnect;
      U_ohm_interconnect = R_interconnect*I_fc;
      U_ohm = U_ohm_electrolyte + U_ohm_interconnect;
//Concentration loss
      U_conc = R*T_cell/(4*F)*(log(p_H2_ch^2*p_O2_ch/(p_H2O_ch^2)) - log(p_H2_eff^2*p_O2_eff/(p_H2O_eff^2)));
      if I_fc>0 then
        R_conc = U_conc/I_fc;
      else
        R_conc = 0;
      end if;
//Activation loss
      i_0 = k_a1*T_cell*exp(-k_a2/(8.314*T_cell));
      U_act1 =  R*T_cell/F*Modelica.Math.asinh(I_fc/ (2*i_0));
      U_act0 = xi_0 + xi_1*T_cell;
      U_act = U_act1 + U_act0;
      if I_fc>0 then
        R_act = U_act1/I_fc;
      else
        R_act = 0;
      end if;
//Double layer charging effect
      U_dl = U_act1 + U_conc;    // U_dl_start;  // U_dl = (I_fc - C_eq*der(U_dl))*(R_act + R_conc);
//Reversible voltage
      E_0 = 1.18 - (T_cell - 298)*k_e;
      E_cell = E_0 + T_cell*R/(4*F)*log(p_H2_ch^2*p_O2_ch/p_H2O_ch^2);
//Output voltage
      U_cell = E_cell - U_act0 - U_ohm - U_dl;
//Thermal submodel
//Air supply tube
//Radiation
      Q_rad_ast = (T_cell^4 - T_ast^4)*A_ast_out*sigma/eps_m_ast_eff;
      eps_m_ast = 1.0807e-10*T_ast^3 - 2.4164e-7*T_ast^2 - 2.2046e-4*T_ast + 0.8962;
      eps_m_ast_eff = 1/eps_m_ast + 87.022/185.34*(1/0.9 -1);
//Outer convection
      Q_out_conv_ast = (T_air_cell - T_ast)*A_ast_out*alpha_ast_out;
      alpha_ast_out = (NusseltNumber_ast_out*k_f_air_out/D_h_ast_out);
      k_f_air_out = 1.5207e-11*T_ast^3 - 4.8574e-8*T_ast^2 + 1.0184e-4*T_ast - 3.9333e-4;
//Inner convection
      Q_in_conv_ast = (T_ast - T_air_ast)*A_ast_in*alpha_ast_in;
      alpha_ast_in = (NusseltNumber_ast_in*k_f_air_in/D_h_ast_in);
      k_f_air_in = 1.5207e-11*T_ast^3 - 4.8574e-8*T_ast^2 + 1.0184e-4*T_ast - 3.9333e-4;
//Heat rate air flow air supply tube
      Q_air_ast = 2*(T_air_ast - T_air_in)*m_air*c_p_air_ast; //T_air_ast_out-T_air_in = 2*(T_air_ast-T_air_in) (s. Simulink)
      m_air = n_c*MM_air;
      c_p_air_ast = 1.9327e-10*T_air_ast^4 - 7.9999e-7*T_air_ast^3 + 1.1407e-3*T_air_ast^2 - 4.489e-1*T_air_ast + 1.0575e3;
//Temperature air in air supply tube
      der(T_air_ast) = Q_in_conv_ast/((c_p_air_ast*(p_c*V_ast_cave/(R*T_air_ast))*MM_air));
//Energy  balance air supply tube
      0 = Q_rad_ast + Q_out_conv_ast - Q_in_conv_ast - Q_air_ast - Q_ast;
//Temperature air supply tube
      der(T_ast) = Q_ast/(c_p_ast*M_ast);
      c_p_ast = 4184*(9.29255e-20*T_ast^6 + 6.58097e-17*T_ast^5 - 1.55425e-12*T_ast^4 + 3.66895e-9*T_ast^3 - 3.79619e-6*T_ast^2 + 1.96313e-3*T_ast - 1.47935e-1);
//Cell
//Heat rate fuel flow
      Q_fuel = (T_fuel_out - T_fuel_in)*C_th_fuel_tot;
      C_th_fuel_tot = C_th_H2_tot + C_th_H2O_tot;
      c_p_H2 = 4184*(3.56903 -4.8950e-4*T_fuel_out + 6.22549e-7*T_fuel_out^2 - 1.19686e-10*T_fuel_out^3) -4124;
      n_H2_tot = n_H2_in + n_H2_out;
      m_H2_tot = n_H2_tot*MM_H2;
      C_th_H2_tot = c_p_H2*m_H2_tot;
      c_p_H2O =  4184*(0.378278 + 1.53443e-4*T_fuel_out + 3.31531e-8*T_fuel_out^2 - 1.78435e-11*T_fuel_out^3)-461.5;
      n_H2O_tot = n_H2O_in + n_H2O_out;
      m_H2O_tot = n_H2O_tot*MM_H2O;
      C_th_H2O_tot = c_p_H2O*m_H2O_tot;
//Convection cell to fuel
      Q_conv_cell_fuel = (T_cell - T_fuel_out)*A_cell*alpha_cell_fuel;
      alpha_cell_fuel = lambda_fuel*NusseltNumber_cell_fuel/D_h_cell_fuel;
      lambda_fuel = lambda_H2 + lambda_H2O;
      K_f_H2 = 0.1*(2.91+(T_fuel_out-600)*0.00345);
      K_f_H2O = 0.1*(0.68+(T_fuel_out-800)*0.001);
      lambda_H2 = K_f_H2*(x_H2_in - (I_fc/(4*F))/n_H2_in);
      lambda_H2O = K_f_H2O*((I_fc/(4*F))/n_H2_in - x_H2_in + 1);
//Heat rate air flow in cell
      Q_air_cell = m_air*c_p_air_cell*(T_air_cell - T_air_ast);
      c_p_air_cell =  1.9327e-10*T_air_cell^4 - 7.9999e-7*T_air_cell^3 + 1.1407e-3*T_air_cell^2 - 4.489e-1*T_air_cell + 1.0575e3;
//Convection cell to air
      Q_conv_cell_air = alpha_cell_air*A_cell*(T_cell - T_air_cell);
      alpha_cell_air = NusseltNumber_cell_air*lambda_air_cell/D_h_cell_air;
      lambda_air_cell = 1.5207e-11*T_cell^3 - 4.8574e-8*T_cell^2 + 1.0184e-4*T_cell - 3.9333e-4;
//Tempterature air in cell
      der(T_air_cell) = Q_conv_cell_air/C_th_air_ann;
      c_p_air_ann = 1.9327e-10*T_air_cell^4 - 7.9999e-7*T_air_cell^3 + 1.1407e-3*T_air_cell^2 - 4.489e-1*T_air_cell + 1.0575e3;
      N_air_ann = V_ann*p_c/(R*T_air_cell);
      M_air_ann = N_air_ann*MM_air;
      C_th_air_ann = M_air_ann*c_p_air_ann;
//Temperature fuel in cell
      der(T_fuel_out) = Q_conv_cell_fuel/(C_th_fuel);
      C_th_fuel = C_th_H2 + C_th_H2O;
      C_th_H2 = c_p_H2*M_H2;
      M_H2 = N_H2*MM_H2;
      N_H2 = p_H2_ch*V_a/(R*T_fuel_out);
      C_th_H2O = c_p_H2O*M_H2O;
      M_H2O = N_H2O*MM_H2O;
      N_H2O = p_H2O_ch*V_a/(R*T_fuel_out);
//Theroetical maximum chemical reaction heat
      Q_max_reac = n_H2_reac*deltah_H2;
//deltah_H2 = 1000*(5.8203e-013*T_cell^4 - 2.2337e-009*T_cell^3 - 1.0957e-007*T_cell^2 + 0.011047*T_cell + 238.75);
      deltah_H2 = deltah_H2_0 + delta_deltah_H2;
//Thermal mass cell
      C_th_cell = c_p_cell*M_cell;
//Electrical power output fuel cell
      P_cell = U_cell*I_fc;
//Energy balance cell
      0 = Q_max_reac - P_cell - Q_rad_ast - Q_conv_cell_air - Q_air_cell - Q_conv_cell_fuel - Q_fuel - Q_cell;
//Cell temperature
      der(T_cell) = Q_cell/C_th_cell;
//Stack variables
      T_stack = T_cell;
      U_stack = U_cell*N_cells;
      E_stack = E_cell*N_cells;
      P_stack = U_stack * I_fc;
      n_H2_in_stack = n_H2_in*N_cells;
      n_H2O_in_stack = n_H2O_in*N_cells;
      n_O2_in_stack = n_O2_in*N_cells;
      m_H2_in_stack = m_H2_in*N_cells;
      m_H2O_in_stack = m_H2O_in*N_cells;
      m_O2_in_stack = m_O2_in*N_cells;
      n_H2_out_stack = n_H2_out*N_cells;
      n_H2O_out_stack = n_H2O_out*N_cells;
      n_O2_out_stack = n_O2_out*N_cells;
      m_H2_out_stack = m_H2_out*N_cells;
      m_H2O_out_stack = m_H2O_out*N_cells;
      m_O2_out_stack = m_O2_out*N_cells;
      if P_stack>0 and n_H2_in>n_H2_out then
         eta_far = (U_cell)/(E_cell)*100;
         eta_el = P_cell/(Q_max_reac)*100;
      else
         eta_far = 0;
         eta_el = 0;
      end if;
//Start variables for integration
//der(I_init) = 0;
//  I_init = 0;
//  p_H2_start =  p_H2_in - I_init*p_a/(4*n_a*F);
//  p_H2O_start =  p_H2O_in + I_init*p_a/(4*n_a*F);
//  p_O2_start = p_O2_in - I_init*p_c/(8*n_c*F);
    annotation(
        experiment(StartTime = 0, StopTime = 9000, Tolerance = 1e-06, Interval = 0.002));
end SOFC;
  equation

  end Components;

  model Records
    extends Modelica.Icons.RecordsPackage;

    partial record PartialRecord
      extends Modelica.Icons.Record;
     
     
     import Modelica.Units.SI;
    //Fuel cell type
    parameter String Type="" annotation (Dialog(group="Type"));
    //Fuel cell specifications
    parameter Integer N_cells                                                       "Number of cells" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Volume V_a                                                         "Volume anode gas channel [m^3]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Volume V_c                                                         "Volume anode gas channel [m^3]" annotation(Dialog(group="Technical Parameters"));
    parameter SI.Length delta                                                       "Length cell [m]" annotation(Dialog(group="Technical Parameters"));
    parameter SI.Area  A_electrolyte                                                "Electrolyte area [m^2]" annotation(Dialog(group="Technical Parameters"));
    parameter SI.Area A_interconnect                                                "Interconnect area [m^2]" annotation(Dialog(group="Technical Parameters"));
    parameter SI.Current I_limit                                                    "Limiting Current [A]" annotation(Dialog(group="Technical Parameters"));
    parameter SI.Temperature T_op                                                   "FC operatingtTemperautre [K]" annotation (Dialog(group="Technical Parameters"));
    
    parameter SI.Area A_ast_out                                                     "Outer area air supply tube [m^2]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Area A_ast_in                                                      "Inner area air supply tube [m^2]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Volume V_ast_cave                                                  "Volume air supply channel cave [m^3]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Mass M_ast                                                         "Mass air supply tube [kg]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Area A_cell                                                        "Area cell [m^2]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Volume V_ann                                                       "Volume of annular space between ast and cell [m^3]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.Mass M_cell                                                        "Mass single cell[kg]" annotation (Dialog(group="Technical Parameters"));
    parameter SI.SpecificHeatCapacity c_p_cell                                      "Specific heat capacity single cell[J/(kgK)]" annotation (Dialog(group="Technical Parameters"));
    //Tunable parameters
    parameter SI.Voltage xi_0                                                       "Tunable parameter [V]" annotation (Dialog(group="Tunable variables"));
    parameter Real xi_1                                                             "Tunable parameter [V/K]" annotation (Dialog(group="Tunable variables"));
    parameter Real k_e                                                              "Tunable parameter [V/K]" annotation (Dialog(group="Tunable variables"));
    parameter Real a_electrolyte                                                    "Tunable parameter" annotation (Dialog(group="Tunable variables"));
    parameter Real b_electrolyte                                                    "Tunable parameter" annotation (Dialog(group="Tunable variables"));
    parameter Real a_interconnect                                                   "Tunable parameter" annotation (Dialog(group="Tunable variables"));
    parameter Real b_interconnect                                                   "Tunable parameter" annotation (Dialog(group="Tunable variables"));
    parameter Real k_a1                                                             "Tunable parameter [A/K]" annotation (Dialog(group="Tunable variables"));
    parameter Real k_a2                                                             "Tunable parameter" annotation (Dialog(group="Tunable variables"));
    //Operating parameters
    parameter SI.Pressure p_a                                                       "Anode pressure [Pa]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.Pressure p_c                                                       "Cathode pressure [Pa]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.MolarFlowRate n_H2_in                                              "Inlet mass flow hydrogen [mol/s]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.MolarFlowRate n_H2O_in                                             "Inlet mass flow water [mol/s]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.MolarFlowRate n_a                                                  "Anode mass flow (fuel flow) [mol/s]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.MolarFlowRate n_c                                                  "Cathode flow (Air flow) [mol/s]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.Temperature T_air_in                                               "Temperature inlet air flow [K]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.Temperature T_fuel_in                                              "Temperature inlet fuel flow [K]" annotation (Dialog(group="Operating Parameters"));
    parameter SI.Temperature T_amb                                                  "Ambient temperature [K]" annotation (Dialog(group="Operating Parameters"));
    parameter Real u_fuel                                                           "Fuel utilization rate as a decimal number"annotation (Dialog(group="Operating Parameters"));
    
      
      
      
    end PartialRecord;

    record Stack_5kW
    
    extends 
        SolidOxideFuelCell.Records.PartialRecord(
        Type="Tubular",
        N_cells = 96,
        V_a = 6.172e-5,
        V_c = 9.902e-5,
        delta = 40e-6,
        A_electrolyte = 200e-4,
        A_interconnect = 45e-4,
        I_limit = 160,
        T_op = 900.2+273.15,
        A_ast_out = 87.022e-4,
        A_ast_in = 62.2035e-4,
        V_ast_cave = 6.1581e-6,
        M_ast = 0.0229882,
        A_cell = 185.34e-4,
        V_ann = 42.6269e-6,
        M_cell = 0.73813,
        c_p_cell = 740,
        xi_0 = 0.15,
        xi_1 = -2e-5,
        k_e = 3.7818e-004,
        a_electrolyte = 2.94e-5,
        b_electrolyte = 10350,
        a_interconnect = 1.256e-3,
        b_interconnect = 4690,
        k_a1 = 53.5041,
        k_a2 = 65000,
        p_a = 3.03975e5,
        p_c = 3.03975e5,
        n_H2_in = 9e-4,
        n_H2O_in = 1e-4,
        n_a = n_H2_in + n_H2O_in,
        n_c = 0.012,
        T_air_in = 1173,
        T_fuel_in = 1173,
        T_amb = 25+273.15,
        u_fuel = 0.85);
    
    

    end Stack_5kW;
  equation

  end Records;
equation

end SolidOxideFuelCell;
