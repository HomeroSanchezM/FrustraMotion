# Charger la molécule
mol new /home/homero/Documentos/M1/S2/Stage/FrustraMotion/results/FRUSTRATION_MTENC/MTENC_CAPSIDS/FRUSTRATION_frames_for_a_monomer/MtEnc_monomer_1_monomers/MtEnc0_monomer1.pdb type pdb waitfor all

# Supprimer les représentations existantes
mol delrep 0 top

# Représentation de base pour toute la molécule
mol representation NewCartoon
mol color Name
mol selection "all"
mol material Opaque
mol addrep top

# Définir les groupes de résidus
set group1 {4 11 19 36 37 39 57 58 69 75 76 78 80 82 84 89 104 111 119 120 128 131 140 149 152 153 171 172 181 198 199 203 204 209 213 214 215 218 222 224 226 230 232 242 244 246 250 252 254 261 263 }
set group2 {10 42 116 129 229 248 }
set group3 {1 2 3 5 6 9 12 13 14 15 16 17 18 20 22 23 24 27 28 29 30 33 34 35 40 41 43 44 45 47 48 49 51 52 53 60 61 62 63 68 71 72 77 81 83 90 95 96 98 100 102 105 106 107 110 115 121 122 124 126 130 133 134 135 136 137 138 144 145 146 151 154 155 159 161 163 167 169 174 175 176 178 179 180 182 184 185 186 192 201 202 205 206 207 208 210 211 212 216 217 219 220 227 228 231 233 235 236 237 239 241 243 247 249 251 255 256 257 259 264 }
set group4 {7 25 26 38 54 55 65 67 70 73 74 79 86 93 99 101 103 108 109 113 118 123 132 139 141 143 148 150 156 158 164 168 170 187 189 190 193 194 196 197 221 223 225 245 260 262 }
set group5 {8 21 31 32 46 50 56 59 64 66 85 87 88 91 92 94 97 112 114 117 125 127 142 147 157 160 162 165 166 173 177 183 188 191 195 200 234 238 240 253 258 }

# Fonction pour créer une sélection et une représentation
proc add_residue_group {resid_list color_id} {
    set selection_text ""
    foreach r $resid_list {
        append selection_text "resid $r or "
    }
    set selection_text [string range $selection_text 0 end-4]

    mol representation NewCartoon
    mol selection $selection_text
    mol color ColorID $color_id
    mol material Opaque
    mol addrep top
}

# apply colors for each group
add_residue_group $group1 19  ;# green2
add_residue_group $group2 1   ;# red
add_residue_group $group3 6   ;# silver
add_residue_group $group4 4   ;# yellow
add_residue_group $group5 5   ;# tan
