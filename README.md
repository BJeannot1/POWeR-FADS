# Program for Observation Well Representation in Fractured Aquifer Dual-continuum Simulations
This repository hosts the source code and example test cases for the Program for Observation Well Representation in Fractured Aquifer Dual-continuum Simulations (POWeR-FADS). This program, intended for use in the context of porous fractured aquifers, aims at simulating the water level in observation wells from fracture and matrix water heads outputed by dual-continuum hydrogeological models. 

It is described more thoroughly in Jeannot et al.(submitted).
If you happen to experience bugs or difficulties in running the code, report it to bjeannot.pro@gmail.com

# I.Introduction
POWeR-FADS postprocesses time series of total water head respectively in the fracture and in the matrix media at the vicinitude of an observation well, as simulated by any planar bidimensional dual-continuum hydrogoelogical model, to simulate the water level in the observation well. POWeR-FADS performs its calculations in a physically based way, by calculating exchange fluxes between the well and each of the two media as a function of head gradients, and introduces as a parameter the altitude of lowest interception of the fracture network by the well. This is a low-parameterized way to describe in the simulations the level of connection between the well and the fractures.
The model is expected to enable researchers better interpret observation wells data of fractured porous aquifers when using bidimensional dual-media hydrogeological models, as it provides a readily available physics-based alternative to the common practice of simply assuming that the well is representative only of fracture water head. 

# **II.	Conceptual model**
![fig1](https://user-images.githubusercontent.com/67539849/194730193-9fae40c2-eaa4-4e1c-92e3-ce047f77d1d0.png)

Fig 1. Conceptual model used in POWeR-FADS for an observation well in a fractured aquifer. Rin,   Rout  and Rdrill are respectively the inner, outer, and drill radius of the observation well. The yellow dots represent the soil of porosity ωdrill used to fill the space between Rout  and Rdrill . ztube and zsurf are respectively the altitude of the top of the well and of the surface. The red plain lines represent the fracture network. Among those, the bold lines represent the saturated portion of the fracture network. The brown and blue backgrounds represent respectively the saturated and unsaturated zones of the matrix. zm and zf are the water heads respectively in the matrix and in the fracture network, and zw is the water level in the well. zbot_m and zbot_f are respectively the altitude of the bottom of the well and the altitude of the lowest interception of the well by the fracture network. The absence of any direct water exchanges between the fracture network and the well below zbot_f has a significant impact on zw. For example, for the represented situation where both zm and zf are below zbot_f, assuming zm and zf are constant in time, then at steady state zw is equal to zm.

For more details on the conceptual model and its undelrying hypotheses, and for explanations about its implementation in the numerical model underpinning POWeR-FADS, see Jeannot et al(submitted).

# **III.	How to use POWeR-FADS**
  ## 1. Requirements

  System : POWeR-FADS uses python and works both on Windows and Linux systems. It has not been tested on Mac.
  
  Python libraries : numpy, matplotlib

  ## 2. Step by step procedure
  - Download the repository
  - Open a command line prompt in the folder containing "POWeR-FADS.py" and type : 
    python POWeR-FADS.py testcase
    
    In the above command, "testcase" is the name of a folder containing an "Inputs" subfolder. The "Inputs" subfolder contains input files corresponding to the test case, as described below. An "Outptus" subfolder will be created in the "testcase" folder as a result of running POWeR-FADS, and outputs of POWeR-FADS will be saved in this subfolder.
    
  For example : "python POWeR-FADS.py Test_case_1" and "python POWeR-FADS.py Test_case_2"  simulate respectively the synthetic test cases 1 and 2 presented in Jeannot et al.(submitted). It is also possible to deal with several test cases in only one call to POWeR-FADS, by writing several command line arguments in sequence when calling the program. Example : "python POWeR-FADS.py Test_case_1 Test_case_2"
    
  ## 3. Description of input files
### zm_zf.txt
The first column is the time in seconds, the second column is total water head in the matrix (zm) in meters, and the third column is the total water head in the fracture network (zf) in meters, as simulated by a dual-continuum bidimensional hydrogeological model at the vicinitude of an observation well. The time step must be constant.
### settings.txt
- write_fluxes : acceptables values are 1 and 0. Output files Fluxes.txt and Sum_fluxes.txt are written only if write_fluxes=1
- write_graphical_outputs : acceptables values are 1 and 0. Output file Output_graph.png is written only if write_graphical_outputs=1
- nbr_Interp_deltaT : The time step at which zm and zf are given might be coarser than the time step at which it is relevant to compute the water level zw in the observation well. nbr_Interp_deltaT is the number of subdivisions to make on the time step of zm and zf, for calculating the time step at which zw is calculated with POWer-FADS. Note that if convergence issues arise, POWeR-FADS will try increasing nbr_Interp_deltaT progressively.
- nbr_Interp_deltaT_MAX : the maximum possible value for nbr_Interp_deltaT
- nbr_pts_rect_int_K : Number of rectangles to use for rectangular integration of the hydraulic conductivity in the vadose zone
- niter_max : maximum number of iterations allowed for convergence at a specific time-step.
- n_failmax : total number of times a set of niter_max iterations can be tried for reaching convergence at a specific-time step. Because POWeR-FADS works with a random relaxation factor from a convergence iteration to the next, convergence can be failed after niter_max iterations on the first try but not on the second one. After n_failmax attemps, the whole POWeR-FADS simulation is restarted from scratch, with an increased nbr_Interp_deltaT.
- epsilon : convergence critetion, in meters.
- Exchanges_well_matrix_through_well_bottom : acceptables values are 1 and 0. If it is set to 1, exchanges between the well and the matrix through the bottom of the well are taken into account in the calculations.
- AveragingK : criterion stating if an arithmetic mean or geometric mean is to use when averaging hydraulic conductivities at the interface between the well and the media, in the case where the well exfiltrates water into the vadose zone. (0=arithmetic,1=geometric)
- calc_auto_zw0 : critetion stating if the initial value of the water level in the well zw should be user-defined, or set automatically by POWeR-FADS from initial values of zf and zm(0=user-defined,1=automatic)
- value_zw0 : user-defined initial value of zw. If calc_auto_zw0=1, this argument is not taken into account by POWeR-FADS
### params_well.txt
- ztube(m): altitude of the top of the tube of the well
- zsurf(m): altitude of the surface
- zbot_f(m) : altitude of lowest interception of the fracture network by the well
- zbot_m(m) : altitude of lthe bottom of the well
- rin(m) : inner radius of the tube
- rout(m) : outer radius of the tube
- rdrill(m) : drill radius of the well
- DeltaE(coupling_length)(m) : width of the interface between the well and the media
### params_matrix.txt
- Ks(m/s) : saturated hydraulic conductivity of the matrix
- n(-) : n coefficient in the Mualem-van Genuchten formula, for the matrix
- alpha(/m): alpha coefficient in the Mualem-van Genuchten formula, for the matrix
- K_unsat_mini(m/s) : lower limit on the saturated hydraulic conductivity calculated with the van Genuchten formula, for the matrix
### params_fractures.txt
- Ks(m/s) : saturated hydraulic conductivity of the fractures
- n(-) : n coefficient in the Mualem-van Genuchten formula, for the fractures
- alpha(/m): alpha coefficient in the Mualem-van Genuchten formula, for the fractures
- K_unsat_mini(m/s) : lower limit on the saturated hydraulic conductivity calculated with the van Genuchten formula, for the fractures

For more details about input formatting, please refer to the example input files provided in the folder "Test_case_1/Inputs" and "Test_case_2/Inputs"

  ## 4. Description of output files
### zw.txt
The first column is the time in seconds and the second column is the water level in the observation well (zw) in meters, as simulated by POWeR-FADS.
### Fluxes.txt
The first column is the time in seconds. The second and third columns are respectively the flux from the well to the matrix and the flux from the well to the fractures in m3/s, as simulated by POWeR-FADS.
### Sum_fluxes.txt
The first column is the time in seconds. The second and third columns are respectively the cumulative flux from the well to the matrix and the cumulative the flux from the well to the fractures in m3, as simulated by POWeR-FADS.
### Output_graph.png
Chart displaying two panels (a) and (b).
Panel (a): Variation over time of total water head in the matrix, noted zm (black dotted line), of total water head in the fracture, noted zf (red dotted line), and of water level in the observation well as simulated by POWeR-FADS, noted zw (blue plain line) . 
Panel (b) : Infiltrated flux from the fracture network to the well (in red) and exfiltrated flux from the well to the matrix (in black), as simulated by POWeR-FADS.

Example for test case 1 :

![Output_graph](https://user-images.githubusercontent.com/67539849/194750310-a5e8e68d-8516-4a3b-94db-c9531fb7f190.png)

# Disclaimer
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
