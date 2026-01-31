# CiliaBuilder CSV Formats

This repository documents the CSV file formats and drawing rules used by **CiliaBuilder** for generating 2D and 3D cilia geometries, including axoneme doublets, central pair microtubules, membrane, cap, cartwheel, and ciliary tip structures.

---

## CiliaBuilder 3D CSV Format

### File format: `3D cilia.csv`

```csv
DoubletNumber,X,Y,Z,Idx_A,Idx_B,Idx_C,Angle,A_Shift,B_Shift,C_Shift
````

### Column definitions

* **DoubletNumber**
  Identifies the structural element represented by the row.

* **X, Y, Z**
  Cartesian coordinates in 3D space.

* **Idx_A, Idx_B**
  Flags indicating whether the A- and/or B-tubule should be drawn.

* **Angle**
  Rotation angle of the tubules.

  * Neighboring doublets differ by 40°
  * Doublet 1 (DMT1) always starts at 90°

* **A_Shift, B_Shift**
  Offset (in Å) from the tubule center to generate the doublet geometry.

  * The same values are used for central pair microtubules

* **C_Shift**
  Structure-specific parameter (see cartwheel section)

---

### DoubletNumber conventions

| DoubletNumber | Structure            |
| ------------: | -------------------- |
|         `1–9` | Doublet microtubules |
|          `-1` | Central pair (CP)    |
|          `-2` | Cap                  |
|          `-3` | Cartwheel            |
|           `0` | Membrane             |

---

### Structure-specific rules

#### Central Pair (CP), Cap, and Membrane

* `Idx_A = 1`, `Idx_B = 1`
* `Angle = 0`

**Central pair shifts**

* `A_Shift = -160`
* `B_Shift = 160`
  (signs may be reversed)

#### Doublet Microtubules

* `A_Shift = -70`
* `B_Shift = 70`

#### Cartwheel (`DoubletNumber = -3`)

* **Angle** → number of cartwheel spokes
* **C_Shift** → centriole radius

---

## CiliaBuilder 2D CSV Format

### Simple centerline format

This format defines only the centerline in 2D. During 3D rendering:

* `Z = Y`
* `Y = 0`

```csv
X,Y
0.0,0.0
…
0.0,100.0
```

---

## Drawing of the Cilia Tip

### Parameters

```text
INITIAL_LENGTH, TRANSITION_LENGTH, TIP_LENGTH
CILIA_RADIUS, TRANSITION_RADIUS, FINAL_RADIUS
```

### Rules

* `INITIAL_LENGTH` enables continuous drawing from the main axoneme

* **B-tubule termination**

  * Random between `INITIAL_LENGTH` and `TRANSITION_LENGTH`
  * Random between `TRANSITION_LENGTH` and `TIP_LENGTH`

* A tip line is generated using a cosine drop:

  * From `INITIAL_LENGTH` to `TIP_LENGTH`
  * Radius decreases from `CILIA_RADIUS` to `FINAL_RADIUS`

* To introduce stochasticity:

  * Generate 9 distinct tip lengths
  * Randomly mix them in memory

---

## Drawing of Primary Cilia

* A template file, `primarycilia_template.csv`, contains all microtubules
* Z coordinates span from `0` to `10000`

### Per-drawing tubule lengths

For each doublet, specify fractional lengths for A- and B-tubules.

**Example**

* `DoubletNumber 1`

  * A-tubule: `0.9`
  * B-tubule: `0.3`

This results in:

* A-tubule drawn to 90% of the total Z length
* B-tubule drawn to 30% of the total Z length

---

### Random drawing rules

* Most B-tubules have length **0.25**
* One B-tubule (e.g. `DoubletNumber 3`) has length **0.40**
* The doublet with `B = 0.40` must have `A = 0.40`
* A-tubule lengths:

  * 4 microtubules: `0.9 – 1.0`
  * 4 microtubules: `0.55 – 0.8`
