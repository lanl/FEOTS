########
Workflow
########

==========================
Overview
==========================
This section of the documentation discusses the steps required to run offline tracer simulations. In short, the following tasks must be accomplished :

* Create impulse fields and passive tracer initial conditions for the parent model
* Run the parent model to produce impulse response functions
* Diagnose transport operators from the impulse and impulse response functions
* Use the transport operators to run offline tracer simulations

For offline tracer simulations, FEOTS has the capability to run in either a global or regional configuration. In the regional configuration, FEOTS provides tools to extract regional transport operators from global operators.

Below, we provide the details for all of the steps in this process while describing how FEOTS global, regional, and simulation databases are created.

Create Transport Operators
==========================
When we discuss "transport operator diagnosis", we specifically mean that we are diagnosing the rows and columns of a sparse matrix that is equivalent to the application of a discrete advection-diffusion operator as implemented in an ocean general circulation model (the parent model). In every time-step of the parent model, we pass a set of impulse fields through passive-tracer forward stepping. The difference in the forward-stepped field and the impulse field, divided by the parent model time step size, provides us with the columns of the transport operator matrix. 

We have provided tools in FEOTS to aid in the process of impulse field generation and operator diagnosis based on these ideas. This allows you to create your own transport operator databases from an online POP model simulation.


Mesh Extraction from a Parent Model
***********************************
.. image:: images/feots_popmesh.png
   :width: 50%
   :alt: FEOTS Mesh Extraction from a parent model

Impulse Field Generation
************************
.. image:: images/feots_impulse.png
   :width: 50%
   :alt: FEOTS Impulse field generation

An impulse field is simply a discrete field of 1's and 0's. Where we place the non-zero impulse values depends on the domain-of-influence of the discrete advection-diffusion operators and the shape of the discrete ocean boundaries. 


Parent Model Execution
***********************

Operator Diagnosis
******************


Run Offline Simulations
=======================

Regional Simulations
====================

Mask Generation
***************
The first step in running a regional simulation is to create a mask file that 

Regional Extraction
*******************

Regional Mapping
****************

Initial Conditions
******************

Forward Integration
*******************

Equilibration
*************


Global Simulations
====================

Initial Conditions
******************

Forward Integration
*******************

Equilibration
*************
