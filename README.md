# generative_design

## Beüzemelés

A környezet teljesen a 3D geometria tárgyban kapott keretre épül. A beüzemelése azzal megegyezik.

## Teszt-modellek
Készítettem Fusion 360-ban egy tesztmodellt, amely `testmodels/` mappában van. Egy adott teszt-modell több különböző fájlból épül fel. Ezért készítettem egy összefogó .txt -t hogy könnyebb legyen a betöltés.
A betöltéshez az elindított tesztprogramon belül a `File/Open` menün belül a load-in.txt kell kiválsztani, amely betölti az összes modellt megfelelően. `(testmodels/2D_01/load-in.txt)`
![open](https://user-images.githubusercontent.com/36598710/94906825-6d993380-049f-11eb-9bf4-a7437aac6255.JPG)

## Használat

A jelenlegi kódban az eredeti környezetben használt gyorsbillentyűk működnek. pl. ` W - wireframed `
Az új funkciókhoz tartozó billentyűk:
- **1** - a kényszer alakzatok megjelenítése/eltűntetése
- **2** - a háromszögek szegmentálása az alapján, hogy melyek tartoznak egyes kényszer alakzatokhoz

![starting_modell](https://user-images.githubusercontent.com/36598710/94906830-6eca6080-049f-11eb-870d-761560f0ae91.JPG)
![afteralgorithm_model](https://user-images.githubusercontent.com/36598710/94906831-6eca6080-049f-11eb-8021-0b617f02d74a.JPG)
