Käyttöesimerkkejä
Suositellaan automaattisen tilan käyttämistä.

Käytä automaattista tilaa. Käsittelee kaikki tiedostot, joiden nimi alkaa avr_force
Lukee voimat force_list.txt nimisestä tiedostosta.
python postprocess.py -a

Käytä manuaalista tilaa. Lue koordinaatteja tiedostoista avr_force0.txt avr_force1.txt ... avr_force19.txt
Lukee voimat force_list.txt nimisestä tiedostosta.
python postprocess.py -n=20

Käytä automaattista tilaa ja tiedostoja, joiden nimet alkavat "voima" esim. voima0.txt
Lukee voimat force_list.txt nimisestä tiedostosta.
python postprocess.py -a -p=voima

Käytä automaattista tilaa. Lue voimat tiedostosta voimat.txt. Tämä on suositeltu vaihtoehto.
python postprocess.py -a -f=voimat.txt

Käytä automaattista tilaa. 50 voimaa generoidaan väliltä [1,10]. Aiheuttaa virhetilanteen,
mikäli askelia on enemmän kuin tiedostoja ja sen takia automaattista tilaa suositellaan.
python postprocess.py -a -f=1,10,50

Käytä manuaalista tilaa. 50 voimaa generoidaan väliltä [1,10]. Luetaan tiedostot, joiden nimet alkavat avr_force. Aiheuttaa virhetilanteen,
mikäli askelia on enemmän kuin tiedostoja ja sen takia automaattista tilaa suositellaan.
python postprocess.py -f=1,10,50

Käytä automaattista tilaa ja voimia, jotka ovat jakautuneet tasaisesti välille [1,10]. Ohjelma laskee voimia välille yhtä paljon kuin on tiedostoja.
python postprocess.py -a -f=1,10
