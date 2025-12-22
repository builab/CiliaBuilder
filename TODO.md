TODO
- Movie part (low priority)
- For ArtiaX build, perhaps presupply a model of doublet & central pair volume at 8-nm repeat (expand curve.py).
- The Doublet now is inverse in order. Might want to change the order but it should not matter much.

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
- Make membrane.
- Group of group modification
- Name change to CiliaBuilder
- Membrane as two layers and cover (Much better) 0.4.4
- When doing sinuisoidal, the CP is actually much shorter. (Priority)
- Make 1 Angstrom different at the starting end as well for better visualization. (Low priority)
- Button to default close the old models & generate new models. 0.4.6
- Change default color
- Make coloring possible for doublet, triplet. 0.4.7
- Make a,b different for cilia as well to make central pair longer 0.4.8
- Generate centriole with offset so it is continuous with cilia
- Make default_config.py for easy controlling of the data. 0.5.2
- Add initial part for sinusoidal to align it. cos(90 - atan2(amplitude, length/(freq*4)))*(2*cilia_radius)
- Use a template curve for line type (csv: X,Y only) and interface update. 0.5.6
- Too many membrane, too many central pair
- generate_tip_curves not ideal yet
- The merged has a wrinkle 0.9.1
- Add to interface & cosine drop better 0.9.3
- Add the cap
- CP rung added
- Add primary to interface, change template to template2D
- Fix Idx_B = 0 for the 1st Z
- Gap in primary cilia
- Change the format of 3D csv to accomodate triplet
- Reorganize read_2D_csv, read_3D_csv to io
- Primary cilia (Using template first)
- Clean up & refactor code (Priority)

