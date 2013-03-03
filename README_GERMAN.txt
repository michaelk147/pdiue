**** Welches Betriebssystem wird unterstützt?

Ubuntu 64 bit 
Linux 3.2.0-27-generic


**** Welche Entwicklungsumgebung brauche ich zum kompilieren der Projekte?

Es wird Eclipse mit CDT verwendet. Eine aktuelle Version sollte funktionieren, ansonsten gilt folgende Version als Referenz:

Eclipse IDE for C/C++ Developers
Version: Juno Release
Build id: 20120614-1722

Als Compiler-Suite kommt GCC zum Einsatz:

gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3


**** Welche Bibliotheken brauche ich?

OpenCV (2.4.2)
Boost (1.50)

gtest (1.60) (nur falls man die Testfälle laufen lassen möchte)


**** Welche Programme brauche ich für das Training?

liblinear-1.91 

**** Wie kompiliere ich die C++ Projekte?

1.) Eclipse mit CDT installieren und starten
2.) Es müssen 2 "path variables" in Eclipse gesetzt werden:
    - Window -> Preferences
    - General -> Workspace -> Linked Resources
    - Dort können über "New..." die Variablen angelegt und mit den lokalen Pfaden belegt werden:
        - LOCAL_BOOST_150  z.B. -> /home/username/builds/boost_1_50_0/install
        - LOCAL_OPENCV_242 z.B. -> /home/username/builds/OpenCV-2.4.2/release/install 
    - Die Struktur dieser Ordner LOCAL_XXX sollte so sein, dass sich dort jeweils direkt "include/" und "lib/" Ordner befinden

3.) Nun können die C++ Projekte aus dem Ordner  "trunk" als Projekte importiert werden.
    Vorgehensweise in Eclipse CDT:
    - File -> Import...
    - General > Existing Projects into Workspace 
    - Click Next
    - Select root directory: Click Browse
	- Hier den Ordner "trunk" auswählen und es werden alle Projekte auch von Unterordnern angezeigt und zum Import angeboten.
    - Checkbox "Copy projects into workspace" setzen
    - Click Finish

    Dieser Import muss für alle folgenden Projekte/Ordner durchgeführt werden:
    - FeatGen
    - LearnDataGenerator
    - pd
    - meanshift/cpp/meanshift  // hier muss der letzte Unterordner im root des "workspace" Ordners liegen.
    - svmrandselect
    - pyramidviz
    
    Falls Fehler von Elclipse angezeigt werden, kann man versuchen über
	 "Rechtsklick auf das Projekt -> Index -> Rebuild"
    den Index neu zu bauen. Dies sollte man direkt einmal für alle Projekte machen (kann dauern).

Builds können nun wie gewohnt über die IDE für jedes Projekte gestartet werden.
    - Project -> Build

Für Release Builds
    - Rechtsklick auf das Projekt -> Build Configurations -> Set Active -> Release
    - Project -> Build
Falls die Includes nicht direkt funktionieren muss Eclipse neugestartet werden.


**** Welche Zusatzprogramme brauche ich für das Training?

liblinear-1.91

Die binaries müssen im path liegen.


**** Welche Zusatzprogramme brauche ich zur Evaluation?

matlab / R2012a


**** Wie starte ich das Training?

Training wird wird über das Shell-Skript "autolearn/autolearn.sh" realisiert. Der Ordner "autolearn" sollte auch auf die lokale Festplatte kopiert werden. Von dort aus kann dann das Skript gestartet werden. In diesem Skript befinden sich weiter Anweisungen/Kommentare.

Kurzfassung:
$bash learnpd.sh rhog6  0.02 2 TESTBUILD

Semantik:
$bash learnpd.sh featid  SVM-C-Parameter NumberOfRetrainings AdditionalNameTag

Damit das Shellscript auf die Release Builds der Programme zugreifen kann, müssen die jeweiligen Projekte in Eclipse als "Release" gebaut werden und entsprechen der Pfad des Eclipses Workspaces in der shell-Datei editiert werden. Für dieses Skript reicht der Build des "LearnDataGenerators" und "svmrandselect" aus.

Nach erfolgreichem Durchlauf befindet sich das trainierte Modell im Ordner "pedestrian_datasets/caltech/data-INRIA/learning/"


**** Wo finde ich bereits trainierte Modelle um den Detector direkt zu starten?

Im Ordner "misc/models" befinden sich alle trainierten Modell-Files. Diese Können dem "pd" Programm übergeben werden ohne dass trainiert werden muss.


**** Wie starte ich den Detektor (pd) ?

Zunächst muss das entsprechende Projekt (pd) in der Release Version über Eclipse gebaut werden.
Danach kann man über die Konsole in den pd Ordner gehen und z.B. über

~/workspaces/cdt/pd$ ./Release/pd -f rhog6 -m ~/pdiue/misc/models/inria.model.rhog6-C-C0.02-B1.RT6 I00000.png -v

Siehe pd Source-Code für weiter Optionen.


**** Wie starte ich die automatische Evaluation ?

Es befindet sich das Shell-Skript "evalmodel.sh" im Ordner "autolearn". Dort finden sich weiter Infos.







