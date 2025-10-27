# occultationFit
Software for fitting 3D asteroid models using stellar occultation observations to determine asteroid sizes.

## Tools <a name="tools"></a>

occultationFit includes several tools described in this section.
The software is modular: if your data is not provided in the required format, you can still build a custom parser that meets the required specifications. The algorithms will still work as expected.

### Event <a name = "event">

This parser handles XML files downloaded from Occult.

REQUIREMENTS:

A custom parser replacing Event must store three dictionaries named star, ast, and observer.

star must be formatted as:

```python
star = {
    "catName": ,
    "catNum": ,
    "gaia_id": ,
    "ra": ,
    "dec":
}
```

- catName: name of the star catalog (e.g., UCAC4)
- catNum: star identifier in the catalog
- gaia_id: Gaia star identifier
- ra, dec: coordinates of the star (right ascension and declination, in degrees)

ast must be formatted as:

```python
ast = {
    "id": ,
    "name":
}
```

- id: asteroid identifier (e.g., 10)
- name: asteroid name (e.g., Hygiea)
Example: (10) Hygiea → id = 10, name = "Hygiea"

observer is a dictionary of dictionaries:
Keys are observer IDs (int), values are:

```python
observer[id] = {
    "name1": ,
    "name2": ,
    "long": ,
    "lat": ,
    "alt": ,
    "isPositive": ,
    "d": Time(),
    "d_err": ,
    "r": Time(),
    "r_err":
}
```

- name1, name2: names of the observer and co-observers
- long, lat, alt: geographical coordinates and elevation of the observing site
- isPositive: boolean indicating whether the chord is positive (occultation detected)
- d, r: disappearance and reappearance timestamps (astropy.time.Time)
- d_err, r_err: timing uncertainties (floats)

A custom parser must also implement:

```python
def __getitem__(self, observer_id):
def __len__(self):
def __iter__(self):
def __next__(self):
```

### ModelParser <a name="model-parser">

This parser handles .tri files and provides multiple methods
(computation of equivalent volume/surface diameter, 3D plotting, inertial tensor calculation and diagonalization, etc.)

REQUIREMENTS:

This parser can be replaced by another one, but the recommended approach is to convert your model into .tri format.

Example of a .tri file:

```text
1018 2032        # number of vertices and number of facets
0.16107606346359549 -6.4287571082150326E-003 1.2520891663492559        # vertex coordinates
0.15112556735743979 4.4585446384259996E-003 1.2500368223958609
0.12954783192959379 2.2386605266567210E-002 1.2461458411015649
...
...
...
3                  # number of vertices per facet
1 2 3              # vertex indices for a facet
3
4 5 6
3
6 5 7
3
6 7 8
```
### Occultation  <a name = "occ">

Take is input an [Event](#event) and manage the construction of the fundamental plan and the placement of cords. It's also offert many functionnalities.


### OccultationFit <a name = "of">

Take as input a list of [Occultation](#occ) and the [ModelParser](#model-parser) for the first pole solution and the second pole solution (if available). This tool is the main tool of this script. Its main functionnality is to scale models according to cords. To do this, the algorithm can use many option from scipy.optimize.minimize as well as a home made solver. It's also possible to combine solvers in way to gain time or to overcome drawbacks and combine advantages.





