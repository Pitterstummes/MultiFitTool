# MultiFitTool

--- Multifittool ---
C Paul Kugler

copy the Multifittool to your datapool folder

Features:
	- No more sorting of files!
	- Automatic detection of all parameters (VTI & dilution fridge)
	- Multifits with Fitmenu
	- Multimodes: Fit multiple modes within a spectra
	- Unique analysis file names with all important information
    - independent of pickorder (first higher, then lower bound works)
	- ...
	- All features of previous fitting tools.

Notes:
    Instead of picking a range with mouseclicking, any key can be used. 
    Just hover the mouse over the desired position and press a key (spacebar,n,..) instead of clicking.

Fitmenu:
	When picking a range for X fits, picking the same point 2 times will open the Fitmenu. (Doubleclick)
	There you can restart the phaserange & fit for that specific sweep (in case of bad fits or other problems),
	end the fitting process for the current sweep (in case of a vanished mode, ...),
	change the multifit number (decreasing multifits as the modes starts to shift) and
	the option to close the menu by pressing X (alt + F4)
	
Speed sometimes (depends on devics) also increases when moving the mouse out of the plot figure.
	
For information regarding possible settings, read comments

UPDATE 12.09.2022:
	- Press (hold down) b while ongoing fit process to get to a break menu in order to repick the range, phase or end the fit.
      If the break menu reopens unwanted, pick ranges with any normal key except 'b' (hover mouse to lower/upper position and
      press button instead of clicking the mouse)
	- NEW SETTING: Fastmode doubles the fitting speed by excluding a figure with real/imaginary parts of the fit.
	
UPDATE 12.10.2022:
	- better structure, new and combined functions, more logic variable names
	- removed irrelevant settings, no more detailed analysis file, no more maxmultifit, no possibility for less process info in figure title
	- new, improved calculation of initial conditions

UPDATE 06.05.2023:
	- gives error if there are wrong "_" in the filenames
	- fixed a error with the current oxford fridge file names
	


	
