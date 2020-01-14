library(colorRamps)
library(RColorBrewer)
library(hash)

# Define modENCODE color spaces:
color.specie.ce = "#66C2A5"
color.specie.dm = "#FC8D62"
color.specie.hs = "#8DA0CB"
color.chromatin = colorRampPalette(c("#f01830", "#f06030", "#f0a830", "#c0c048", "#90c048", "#78c060", "#60c078", "#009048", "#186030", "#784890", "#a84890", "#787878", "#a8a8a8", "#f0c0d8", "#fff0f0", "#f0f0d8"), space="Lab")

# Define a function to reverse color ramps:
revColor <- function(color.palette) { return(function(x)rev(color.palette(x))) }

# hand-made colors, in a rush.
citric = colorRampPalette(c("#03AB11", "#FFF301"), space="Lab")
citrus = colorRampPalette(c("#15E602", "#FFFF33"), space="Lab")
berry = colorRampPalette(c("#7700FF","#FF0080"), space="Lab")
forest = colorRampPalette(c("#E0E0E0","#6BDA61", "#15E602"), space="Lab")
white.mango = colorRampPalette(c("#FFFFFF", "#FF00FF", "#FF0000"), space="Lab")
white.tango = colorRampPalette(c("#FFFFFF", "#FF00FF", "#9500FF"), space="Lab")
white.grove = colorRampPalette(c("#FFFFFF", "#FFF301", "#BAE61D", "#01DF29"), space="Lab")
white.jungle = colorRampPalette(c("#FFFFFF", "#FFF301","#03AB11"), space="Lab")
white.orange = colorRampPalette(c("#FFFFFF", "#FFF301","#FF7700"), space="Lab")
horizon = colorRampPalette(c("#000033", "#000075", "#0000B6", "#0000F8", "#2E00FF", "#6100FF", "#9408F7", "#C729D6", "#FA4AB5", "#FF6A95", "#FF8B74", "#FFAC53", "#FFCD32", "#FFEE11", "#FFFF60"), space="Lab")

# hand-made colors, from R.
wolfgang.basic = colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"), space="Lab")
wolfgang.extra = colorRampPalette(c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", "#2A9EC1", "#206AAD", "#243996", "#081D58"), space="Lab")
solar.flare = colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"), space="Lab")
solar.glare = colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FCFCFC", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"), space="Lab")
solar.basic = colorRampPalette(c("#214B85", "#1873CC", "#1E90FF", "#00BFFF", "#ACD8E5", "#D2D2D2", "#FFD700", "#ED2C2C", "#A31D1D"), space="Lab")
solar.extra = colorRampPalette(c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D"), space="Lab")
solar.blues = colorRampPalette(c("#FCFCFC", "#C0E4FD", "#75CEFE", "#0CB9FF", "#1BA7FF", "#1E95FF", "#2884E7", "#3072C5", "#3361A5"), space="Lab")
solar.rojos = colorRampPalette(c("#FCFCFC", "#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", "#F14C2B", "#DA2828", "#BE2222", "#A31D1D"), space="Lab")
samba.color = colorRampPalette(c("#1B85ED", "#1AA2F6", "#00BFFF", "#4AC596", "#00CC00", "#A7D400", "#FFD700", "#FFBE00", "#FFA500"), space="Lab")
samba.night = colorRampPalette(c("#1873CC", "#1798E5", "#00BFFF", "#4AC596", "#00CC00", "#A2E700", "#FFFF00", "#FFD200", "#FFA500"), space="Lab")
samba.light = colorRampPalette(c("#00D1FF","#03AB11","#FFF301"), space="Lab")

# hand-made colors, from show.
dusk.dawn = colorRampPalette(c("#98ABC5", "#8D91AD", "#827896", "#775F80", "#6B476B", "#93575B", "#B8684A", "#DB7933", "#FF8C00"), space="Lab")
dark.cyan = colorRampPalette(c("#000000", "#0E2824", "#014C44", "#0A5F4F", "#13725A", "#19997F", "#1EC0A6", "#19DFD2", "#00FFFF"), space="Lab")
dark.blue = colorRampPalette(c("#000000", "#00171F", "#002F3F", "#00475F", "#005F7F", "#00779F", "#008FBF", "#00A7DF", "#00BFFF"), space="Lab")
dark.citrus = colorRampPalette(c("#000000", "#22350F", "#3B680C", "#529111", "#6ABB15", "#74DD0F", "#7FFF00", "#ADF121", "#D1E131"), space="Lab")
dark.violet = colorRampPalette(c("#000000", "#1E0A35", "#31016A", "#4B0181", "#660099", "#7800CA", "#8A00FF", "#C800FF", "#FE00FF"), space="Lab")
ocean.green = colorRampPalette(c("#07519B", "#2975B4", "#5097C9", "#93C1DF", "#FCFCFC", "#CAEAC5", "#97D494", "#5BAB5A", "#006400"), space="Lab")
ocean.earth = colorRampPalette(c("#0F3341", "#1563AA", "#0B99E6", "#3DCDFD", "#F7F7F7", "#B87350", "#872E1C", "#601622", "#401C2A"), space="Lab")
ocean.brick = colorRampPalette(c("#0F3341", "#1563AA", "#0B99E6", "#3DCDFD", "#F7F7F7", "#EB9457", "#D1551F", "#B02F1B", "#8D1616"), space="Lab")
algae.earth = colorRampPalette(c("#543005", "#985D12", "#CFA154", "#F0DEB1", "#F5F5F5", "#B5E2DC", "#5AB2A8", "#0E726A", "#003C30"), space="Lab")
flame.flame = colorRampPalette(c("#000033", "#0000A5", "#1E00FB", "#6F00FD", "#C628D6", "#FE629D", "#FF9B64", "#FFD52C", "#FFFF5F"), space="Lab")
flame.light = colorRampPalette(c("#000033", "#000E92", "#1300FF", "#8E0EEA", "#C628D6", "#E9699F", "#FF9B63", "#FFCE62", "#FFFF5F"), space="Lab")
flame.polar = colorRampPalette(c("#C628D6", "#8E0EEA", "#1300FF", "#000E92", "#000033", "#7F494D", "#FF9B63", "#FFCE62", "#FFFF5F"), space="Lab")
flame.volts = colorRampPalette(c("#000000", "#371377", "#5F00FF", "#9400FF", "#BE00FF", "#E000EB", "#FF00D8", "#FF0090", "#FF004B"), space="Lab")
flame.watts = colorRampPalette(c("#FFFFFF", "#C190FF", "#5F00FF", "#9400FF", "#BE00FF", "#E000EB", "#FF00D8", "#FF0090", "#FF004B"), space="Lab")
flame.artic = colorRampPalette(c("#000000", "#371377", "#5F00FF", "#BD00EC", "#FF00D8", "#C7ACEC", "#00FFFF", "#0AD7D3", "#0DB2AA"), space="Lab")
flame.weird = colorRampPalette(c("#00FFFF", "#0AD7D3", "#0DB2AA", "#1C5551", "#000000", "#371377", "#5F00FF", "#BD00EC", "#FF00D8"), space="Lab")
flame.blind = colorRampPalette(c("#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF", "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"), space="Lab")
flame.macaw = colorRampPalette(c("#000000", "#28410F", "#477C0E", "#64B114", "#9FCF23", "#C9E553", "#81F7D0", "#16DCD2", "#1AA58C"), space="Lab")
flame.wings = colorRampPalette(c("#D1E131", "#85C51D", "#529111", "#2F4E0F", "#000000", "#0F4338", "#107E6A", "#1BBBA7", "#00FFFF"), space="Lab")

# following colors are generated online (http.//www.pixelfor.me/crc/).
calma.azules = colorRampPalette(c("#031C25", "#093B4D", "#1C5F77", "#3685A2", "#56A6C3", "#86C2D8", "#B6DDEB", "#F2FBFE"), space="Lab")
calma.musgos = colorRampPalette(c("#212503", "#444D09", "#6B771C", "#93A236", "#B4C356", "#CDD886", "#E4EBB6", "#FCFEF2"), space="Lab")
calma.bosque = colorRampPalette(c("#032506", "#094D0E", "#1C7722", "#36A23D", "#56C35D", "#86D88B", "#B6EBBA", "#F2FEF3"), space="Lab")
calma.marino = colorRampPalette(c("#032515", "#094D2D", "#1C774D", "#36A26F", "#56C390", "#86D8B2", "#B6EBD2", "#F2FEF8"), space="Lab")
calma.morado = colorRampPalette(c("#030925", "#09154D", "#1C2B77", "#3648A2", "#5668C3", "#8694D8", "#B6BFEB", "#F2F4FE"), space="Lab")
calma.manudo = colorRampPalette(c("#290303", "#590707", "#8C1616", "#BE2A2A", "#DF4A4A", "#ED8080", "#F7B4B4", "#FFEEEE"), space="Lab")
china.theory = colorRampPalette(c("#120324", "#420A4A", "#721D57", "#9B3850", "#BC6B58", "#D3B687", "#E6E8B7", "#F8FDF2"), space="Lab")
china.ranges = colorRampPalette(c("#031424", "#1F0A4A", "#721D64", "#9B3838", "#BCAB58", "#A0D387", "#B7E8CF", "#F2F9FD"), space="Lab")
china.weirdo = colorRampPalette(c("#04032E", "#2E0267", "#890BA3", "#DE15AF", "#FF347E", "#FF7772", "#FFCFAB", "#FFFBEA"), space="Lab")
china.basics = colorRampPalette(c("#25032E", "#670253", "#A30B48", "#DE1515", "#FF8534", "#FFE272", "#EEFFAB", "#F2FFEA"), space="Lab")
china.sunset = colorRampPalette(c("#031124", "#0C0A4A", "#451D72", "#91389B", "#BC589B", "#D38799", "#E8C1B7", "#FDF9F2"), space="Lab")
china.dragon = colorRampPalette(c("#03032A", "#2B065C", "#801491", "#C52696", "#E74671", "#F2917D", "#FADEB3", "#FEFFED"), space="Lab")
china.novice = colorRampPalette(c("#2A0E03", "#5C4406", "#7E9114", "#68C526", "#46E748", "#7DF2B2", "#B3FAF2", "#EDF9FF"), space="Lab")

# following colors are from R ColorBrewer.
brewer.fire = colorRampPalette(c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"), space="Lab")
brewer.heat = colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"), space="Lab")
brewer.orange = colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"), space="Lab")
brewer.red = colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"), space="Lab")
brewer.green = colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), space="Lab")
brewer.blue = colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"), space="Lab")
brewer.purple = colorRampPalette(c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"), space="Lab")
brewer.violet = colorRampPalette(c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"), space="Lab")
brewer.jamaica = colorRampPalette(c("#006837", "#2DA154", "#86CB66", "#CCE982", "#FFFFBF", "#FDD380", "#F88D51", "#DE3F2E", "#A50026"), space="Lab")
brewer.marine = colorRampPalette(c("#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"), space="Lab")
brewer.spectra = colorRampPalette(c("#5E4FA2", "#3F96B7", "#88CFA4", "#D7EF9B", "#FFFFBF", "#FDD380", "#F88D51", "#DC494C", "#9E0142"), space="Lab")
brewer.celsius = colorRampPalette(c("#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF", "#FDD384", "#F88D51", "#DE3F2E", "#A50026"), space="Lab")
brewer.yes = colorRampPalette(c("#053061", "#2971B1", "#6AACD0", "#C1DDEB", "#F7F7F7", "#FACDB5", "#E58267", "#BB2933", "#67001F"), space="Lab")

# following colors generated from http.//www.colorhexa.com/examples
forest.yellow = colorRampPalette(c("#215e33", "#306835", "#3f7136", "#4e7b38", "#5d843a", "#6c8e3b", "#7b983d", "#89a13f", "#98ab41", "#a7b542", "#b6be44", "#c5c846", "#d4d147", "#e3db49"), space="Lab")
forest.citric = colorRampPalette(c("#0b3310", "#1b4210", "#2b5111", "#3b6111", "#4c7012", "#5c7f12", "#6c8e12", "#7c9e13", "#8cad13", "#9cbc13", "#adcb14", "#bddb14", "#cdea15", "#ddf915"), space="Lab")
citric.yellow = colorRampPalette(c("#21be45", "#30c246", "#40c546", "#4fc947", "#5fcc48", "#6ed048", "#7dd349", "#8dd74a", "#9cda4b", "#abde4b", "#bbe14c", "#cae54d", "#dae84d", "#e9ec4e"), space="Lab")
ocean.citrus = colorRampPalette(c("#3683ba", "#418bb0", "#4c93a7", "#569b9d", "#61a393", "#6cab8a", "#77b380", "#81ba76", "#8cc26c", "#97ca63", "#a2d259", "#acda4f", "#b7e246", "#c2ea3c"), space="Lab")
ocean.pink = colorRampPalette(c("#3b5e84", "#4a5e84", "#595e84", "#685e84", "#765e84", "#855e84", "#945e84", "#a35f85", "#b25f85", "#c15f85", "#cf5f85", "#de5f85", "#ed5f85", "#fc5f85"), space="Lab")
ocean.red = colorRampPalette(c("#3b82ae", "#487aa1", "#547294", "#616987", "#6d617a", "#7a596d", "#865160", "#934853", "#9f4046", "#ac3839", "#b8302c", "#c5271f", "#d11f12", "#de1705"), space="Lab")
ocean.aqua = colorRampPalette(c("#2668aa", "#2d6eaa", "#3574aa", "#3c7aaa", "#4380a9", "#4b86a9", "#528ca9", "#5993a9", "#6099a9", "#689fa9", "#6fa5a8", "#76aba8", "#7eb1a8", "#85b7a8"), space="Lab")
ocean.teal = colorRampPalette(c("#0a0a66", "#15176c", "#202573", "#2b3279", "#364080", "#414d86", "#4c5a8d", "#576893", "#62759a", "#6d82a0", "#7890a7", "#839dad", "#8eabb4", "#99b8ba"), space="Lab")
cyan.brick = colorRampPalette(c("#6cd0c2", "#70c3b3", "#74b6a5", "#78a996", "#7b9c88", "#7f8f79", "#83826a", "#87755c", "#8b684d", "#8f5b3e", "#924e30", "#964121", "#9a3413", "#9e2704"), space="Lab")
aqua.brick = colorRampPalette(c("#019bcf", "#0d94c2", "#1a8cb4", "#2685a7", "#337d99", "#3f768c", "#4c6e7e", "#586771", "#655f63", "#715856", "#7e5048", "#8a493b", "#97412d", "#a33a20"), space="Lab")
aqua.tan = colorRampPalette(c("#62adbb", "#6cb0b5", "#75b2af", "#7fb5a8", "#88b7a2", "#92ba9c", "#9bbd96", "#a5bf8f", "#aec289", "#b8c583", "#c1c77d", "#cbca76", "#d4cc70", "#decf6a"), space="Lab")
cyan.tan = colorRampPalette(c("#4fe8c2", "#5be3bb", "#67deb5", "#73d9ae", "#7fd3a7", "#8bcea1", "#97c99a", "#a4c493", "#b0bf8c", "#bcba86", "#c8b47f", "#d4af78", "#e0aa72", "#eca56b"), space="Lab")
teal.orange = colorRampPalette(c("#0cb499", "#1dae92", "#2ea98b", "#3fa384", "#4f9e7d", "#609876", "#71936f", "#828d69", "#938862", "#a4825b", "#b47d54", "#c5774d", "#d67246", "#e76c3f"), space="Lab")
teal.violet = colorRampPalette(c("#4b8b84", "#508184", "#567784", "#5b6d84", "#616384", "#665984", "#6b4f84", "#714683", "#763c83", "#7b3283", "#812883", "#861e83", "#8c1483", "#910a83"), space="Lab")
blue.cyan = colorRampPalette(c("#4111f2", "#4923f0", "#5135ed", "#5946eb", "#6258e8", "#6a6ae6", "#727ce3", "#7a8de1", "#829fde", "#8ab1dc", "#93c3d9", "#9bd4d7", "#a3e6d4", "#abf8d2"), space="Lab")
purple.pink = colorRampPalette(c("#6848d1", "#7345ca", "#7f42c3", "#8a3fbc", "#953db5", "#a13aae", "#ac37a7", "#b734a1", "#c2319a", "#ce2e93", "#d92c8c", "#e42985", "#f0267e", "#fb2377"), space="Lab")
purple.baby = colorRampPalette(c("#511293", "#5a239b", "#6234a2", "#6b45aa", "#7356b2", "#7c67b9", "#8478c1", "#8d88c9", "#9599d1", "#9eaad8", "#a6bbe0", "#afcce8", "#b7ddef", "#c0eef7"), space="Lab")
cyan.green = colorRampPalette(c("#29ddea", "#31dcd8", "#39dcc7", "#41dbb5", "#4adba3", "#52da92", "#5ad980", "#62d96e", "#6ad85c", "#72d74b", "#7bd739", "#83d627", "#8bd616", "#93d504"), space="Lab")
cyan.pink = colorRampPalette(c("#606bef", "#6c6ae6", "#7869dd", "#8468d4", "#9168cb", "#9d67c2", "#a966b9", "#b565b1", "#c164a8", "#cd639f", "#da6396", "#e6628d", "#f26184", "#fe607b"), space="Lab")
cyan.violet = colorRampPalette(c("#29e9ae", "#36dbae", "#43cdae", "#50beaf", "#5eb0af", "#6ba2af", "#7894af", "#8585b0", "#9277b0", "#9f69b0", "#ad5bb0", "#ba4cb1", "#c73eb1", "#d430b1"), space="Lab")
cyan.purple = colorRampPalette(c("#4adbf1", "#51cff0", "#59c3ef", "#60b7ee", "#68abed", "#6f9fec", "#7793eb", "#7e88eb", "#867cea", "#8d70e9"), space="Lab")


# Define 1st generation color spaces:
color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")
color.palette1 = colorRampPalette(c("#3362A5", "dodgerblue1", "deepskyblue", "white", "gold", "firebrick2","#A31D1D"), space="Lab")
color.palette2 = colorRampPalette(c("#3362A5", "dodgerblue1", "deepskyblue", "gray99", "gold", "firebrick2","#A31D1D"), space="Lab")
color.palette3 = colorRampPalette(c("gray99", "deepskyblue", "dodgerblue1", "#3362A5"), space="Lab")
color.palette4 = colorRampPalette(c("gray99", "gold", "firebrick2","#A31D1D"), space="Lab")
color.palette5 = colorRampPalette(c("white",brewer.pal(9, "YlGnBu")))
color.palette6 = colorRampPalette(c("#224B85", "dodgerblue3", "dodgerblue1", "deepskyblue", "lightblue", "lightgray", "gold", "firebrick2","#A31D1D"), space="Lab")
color.palette7 = colorRampPalette(c("#3362A5", "dodgerblue1", "deepskyblue", "lightblue", "lightgray", "gold", "firebrick2","#A31D1D"), space="Lab")
color.palette8 = colorRampPalette(c("dodgerblue2", "deepskyblue", "green3", "gold", "orange"), space="Lab")
color.palette9 = colorRampPalette(c("dodgerblue3", "deepskyblue", "green3", "yellow", "orange"), space="Lab")
color.paletteR = colorRampPalette(brewer.pal(9,"Reds"))
color.paletteG = colorRampPalette(brewer.pal(9,"Greens"))
color.paletteB = colorRampPalette(brewer.pal(9,"Blues"))
color.paletteO = colorRampPalette(brewer.pal(9,"Oranges"))
color.paletteP = colorRampPalette(brewer.pal(9,"Purples"))
color.paletteA = colorRampPalette(brewer.pal(9,"BuGn")) #Algae
color.paletteM = colorRampPalette(brewer.pal(9,"GnBu")) #Marine
color.paletteU = colorRampPalette(brewer.pal(9,"BuPu")) #UV
color.paletteD = colorRampPalette(brewer.pal(9,"PuBu")) #Death
color.paletteH = colorRampPalette(brewer.pal(9,"OrRd")) #Heat
color.paletteV = colorRampPalette(brewer.pal(9,"RdPu")) #Violet
color.paletteC = colorRampPalette(brewer.pal(9,"YlGn")) #Citric
color.paletteW = colorRampPalette(brewer.pal(9,"YlGnBu")) #Wolfgang
color.paletteF = colorRampPalette(brewer.pal(9,"YlOrBr")) #Fire
color.paletteN = colorRampPalette(brewer.pal(11,"BrBG")) #Natural
color.paletteE = colorRampPalette(brewer.pal(11,"PiYG")) #Extreme
color.paletteI = colorRampPalette(brewer.pal(11,"PRGn")) #Irradiate
color.paletteZ = colorRampPalette(brewer.pal(11,"PuOr")) #Zigzag
color.paletteY = revColor(colorRampPalette(brewer.pal(11,"RdBu"))) #Yes!
color.paletteL = colorRampPalette(brewer.pal(11,"RdGy")) #Liga
color.paletteT = revColor(colorRampPalette(brewer.pal(11,"RdYlBu"))) #Temperature/Transition
color.paletteJ = revColor(colorRampPalette(brewer.pal(11,"RdYlGn"))) #Jamaica
color.paletteS = revColor(colorRampPalette(brewer.pal(11,"Spectral"))) #Spectral
color.paletteQ = colorRampPalette(c("dodgerblue3", "deepskyblue", "lightgray", "gold","orange", "firebrick2", "firebrick3"), space="Lab")
color.paletteK = colorRampPalette(c("dodgerblue3", "dodgerblue1", "deepskyblue", "lightblue", "lightgray", "gold", "firebrick2","firebrick3"), space="Lab")
color.paletteX = colorRampPalette(c("#3362A5", "dodgerblue3", "dodgerblue2", "dodgerblue1", "deepskyblue", "lightblue", "lightgray", "gold", "orange","firebrick2"), space="Lab")
color.paletteA1 = colorRampPalette(c("black","deepskyblue")) #DeepSkyBlue
color.paletteA2 = colorRampPalette(c("#08519c", "#3182bd", "#6baed6", "gray99", "#bae4b3", "#74c476", "darkgreen"), space="Lab")
color.paletteA3 = colorRampPalette(c("#0f3341", "dodgerblue3", "deepskyblue", "gray97", "#a44819", "#6b1520", "#401c2a"))
color.paletteA4 = colorRampPalette(c("#0f3341", "dodgerblue3", "deepskyblue", "gray97", "#E77422", "#BC381D", "#8D1616"))
color.paletteA5 = colorRampPalette(c("#1C3198", "#8DA0CB", "gray99", "#66C2A5", "#0A7511"), space="Lab")
color.paletteA6 = colorRampPalette(c("#1C3198", "#8DA0CB", "gray99", "#FC8D62", "#A31D1D"), space="Lab")

# List all color palettes:
color.palettes = hash()
color.palettes["color.palette0"]=color.palette0
color.palettes["color.palette1"]=color.palette1
color.palettes["color.palette2"]=color.palette2
color.palettes["color.palette3"]=color.palette3
color.palettes["color.palette4"]=color.palette4 
color.palettes["color.palette5"]=color.palette5
color.palettes["color.palette6"]=color.palette6
color.palettes["color.palette7"]=color.palette7
color.palettes["color.palette8"]=color.palette8
color.palettes["color.palette9"]=color.palette9
color.palettes["color.paletteR"]=color.paletteR
color.palettes["color.paletteG"]=color.paletteG
color.palettes["color.paletteB"]=color.paletteB
color.palettes["color.paletteP"]=color.paletteP
color.palettes["color.paletteO"]=color.paletteO
color.palettes["color.paletteA"]=color.paletteA
color.palettes["color.paletteM"]=color.paletteM
color.palettes["color.paletteU"]=color.paletteU
color.palettes["color.paletteD"]=color.paletteD
color.palettes["color.paletteH"]=color.paletteH
color.palettes["color.paletteV"]=color.paletteV
color.palettes["color.paletteC"]=color.paletteC
color.palettes["color.paletteW"]=color.paletteW
color.palettes["color.paletteF"]=color.paletteF
color.palettes["color.paletteN"]=color.paletteN
color.palettes["color.paletteE"]=color.paletteE
color.palettes["color.paletteI"]=color.paletteI
color.palettes["color.paletteZ"]=color.paletteZ
color.palettes["color.paletteY"]=color.paletteY
color.palettes["color.paletteL"]=color.paletteL
color.palettes["color.paletteT"]=color.paletteT
color.palettes["color.paletteJ"]=color.paletteJ
color.palettes["color.paletteS"]=color.paletteS
color.palettes["color.paletteQ"]=color.paletteQ
color.palettes["color.paletteK"]=color.paletteK
color.palettes["color.paletteX"]=color.paletteX
color.palettes["color.paletteA1"]=color.paletteA1
color.palettes["color.paletteA2"]=color.paletteA2
color.palettes["color.paletteA3"]=color.paletteA3
color.palettes["color.paletteA4"]=color.paletteA4
color.palettes["color.paletteA5"]=color.paletteA5
color.palettes["color.paletteA6"]=color.paletteA6

# Define 2nd generation color spaces:
color.parrotAX = colorRampPalette(c("#000000", "#31026A", "#660099", "#8A00FF","#FE00FF"), space="Lab")
color.parrotAY = colorRampPalette(c("#000033FF", "#000075FF", "#0000B6FF", "#0000F8FF", "#2E00FFFF", "#6100FFFF", "#9408F7FF", "#C729D6FF", "#FA4AB5FF", "#FF6A95FF", "#FF8B74FF", "#FFAC53FF", "#FFCD32FF", "#FFEE11FF", "#FFFF60FF"), space="Lab")
color.parrotA0 = colorRampPalette(c("#000033FF", "#1400FFFF", "#C729D6FF", "#FF9C63FF", "#FFFF60FF"), space="Lab")
color.parrotA1 = colorRampPalette(c("#C729D6FF", "#1400FFFF", "#000033FF", "#FF9C63FF", "#FFFF60FF"), space="Lab")
color.parrotA2 = colorRampPalette(c("#C729D6FF", "#1400FFFF", "#000033FF", "#FF9C63FF", "#FFFF60FF"), space="Lab")
color.parrotA3 = colorRampPalette(c("#FFFFFF", "#1400FFFF", "#C729D6FF", "#FF9C63FF", "#FFFF60FF"), space="Lab")
color.parrotA4 = colorRampPalette(c("#1400FFFF", "#C729D6FF", "#FFFFFF", "#FFFF60FF", "#FF9C63FF"), space="Lab")
color.parrotA5 = colorRampPalette(c("#1400FFFF", "#C729D6FF", "#FFFFFF", "#FFFF60FF", "#FF9C63FF"), space="Lab")

color.palettes["color.parrotAX"]=color.parrotAX
color.palettes["color.parrotAY"]=color.parrotAY
color.palettes["color.parrotA0"]=color.parrotA0
color.palettes["color.parrotA1"]=color.parrotA1
color.palettes["color.parrotA2"]=color.parrotA2
color.palettes["color.parrotA3"]=color.parrotA3
color.palettes["color.parrotA4"]=color.parrotA4
color.palettes["color.parrotA5"]=color.parrotA5
#printColors(c("color.parrotA4"))

color.parrotBX = colorRampPalette(c("#000000", "#3C690C", "#6ABC16", "#7FFF00", "#D1E231"), space="Lab")
color.parrotBY = colorRampPalette(c("#000000", "#024D44", "#14735B", "#1EC1A7", "#00FFFF"), space="Lab")
color.parrotB0 = colorRampPalette(c("#000000", "#3C690C", "#6ABC16", "#D1E231", "#00FFFF", "#1BA68C"), space="Lab")
color.parrotB1 = colorRampPalette(c("#D1E231", "#6ABC16", "#3C690C", "#000000", "#07594B", "#1BA68C", "#00FFFF"), space="Lab")
color.parrotB2 = colorRampPalette(c("#D1E231", "#6ABC16", "#3C690C", "#000000", "#07594B", "#1BA68C", "#00FFFF"), space="Lab")
color.parrotB3 = colorRampPalette(c("#FFFFFF", "#D1E231", "#6ABC16", "#3C690C", "#07594B", "#1BA68C", "#00FFFF"), space="Lab")
color.parrotB4 = colorRampPalette(c("#3C690C", "#6ABC16", "#D1E231", "#FFFFFF", "#00FFFF", "#1BA68C", "#07594B"), space="Lab")
color.parrotB5 = colorRampPalette(c("#3C690C", "#6ABC16", "#D1E231", "#FFFFFF", "#00FFFF", "#1BA68C", "#07594B"), space="Lab")

color.palettes["color.parrotBX"]=color.parrotBX
color.palettes["color.parrotBY"]=color.parrotBY
color.palettes["color.parrotB0"]=color.parrotB0
color.palettes["color.parrotB1"]=color.parrotB1
color.palettes["color.parrotB2"]=color.parrotB2
color.palettes["color.parrotB3"]=color.parrotB3
color.palettes["color.parrotB4"]=color.parrotB4
color.palettes["color.parrotB5"]=color.parrotB5
#printColors(c("color.parrotB0"))

color.parrotC0 = colorRampPalette(c("#000000", "#6000FF", "#BE00FF", "#FF00D9", "#FF004B"), space="Lab")
color.parrotC1 = colorRampPalette(c("#6000FF", "#BE00FF", "#000000", "#FF00D9", "#FF004B"), space="Lab")
color.parrotC2 = colorRampPalette(c("#4F2398", "#6000FF", "#BE00FF", "#000000", "#FF00D9", "#FF004B", "#B20931"), space="Lab")
color.parrotC3 = colorRampPalette(c("#FFFFFF", "#6000FF", "#BE00FF", "#FF00D9", "#FF004B"), space="Lab")
color.parrotC4 = colorRampPalette(c("#6000FF", "#BE00FF", "#FFFFFF", "#FF00D9", "#FF004B"), space="Lab")
color.parrotC5 = colorRampPalette(c("#4F2398", "#6000FF", "#BE00FF", "#FFFFFF", "#FF00D9", "#FF004B", "#B20931"), space="Lab")

color.palettes["color.parrotC0"]=color.parrotC0
color.palettes["color.parrotC1"]=color.parrotC1
color.palettes["color.parrotC2"]=color.parrotC2
color.palettes["color.parrotC3"]=color.parrotC3
color.palettes["color.parrotC4"]=color.parrotC4
color.palettes["color.parrotC5"]=color.parrotC5
#printColors(c("color.parrotC0"))

color.parrotD0 = colorRampPalette(c("#000000", "#6000FF", "#FF00D9", "#00FFFF", "#0EB2AA"), space="Lab")
color.parrotD1 = colorRampPalette(c("#00FFFF", "#0EB2AA", "#000000", "#6000FF", "#FF00D9"), space="Lab")
color.parrotD2 = colorRampPalette(c("#00FFFF", "#0EB2AA", "#000000", "#6000FF", "#FF00D9"), space="Lab")
color.parrotD3 = colorRampPalette(c("#FFFFFF", "#FF00D9", "#6000FF", "#0EB2AA", "#00FFFF"), space="Lab")
color.parrotD4 = colorRampPalette(c("#0EB2AA", "#00FFFF", "#FFFFFF", "#FF00D9", "#6000FF"), space="Lab")
color.parrotD5 = colorRampPalette(c("#0EB2AA", "#00FFFF", "#FFFFFF", "#FF00D9", "#6000FF"), space="Lab")

color.palettes["color.parrotD0"]=color.parrotD0
color.palettes["color.parrotD1"]=color.parrotD1
color.palettes["color.parrotD2"]=color.parrotD2
color.palettes["color.parrotD3"]=color.parrotD3
color.palettes["color.parrotD4"]=color.parrotD4
color.palettes["color.parrotD5"]=color.parrotD5
#printColors(c("color.parrotD5"))


color.parrotE0 = colorRampPalette(c("white", "deepskyblue", "green3", "yellow", "orange"), space="Lab")

reScale <- function(values, mode="fraction", percent=FALSE) {
  if (mode == "fraction") { values = values-min(values) ; values = values/max(values) }
  else if (mode == "max") { values = values/max(values) }
  else if (mode == "min") { values = values/min(values) }
  if (percent) { values = 100*values }
  return(values)
}

convertZeros <- function(vector, min.value="minimum") {
  if (min.value == "minimum") { min.value = min(vector[vector !=0]) }
  output = sub(0, min.value, vector)
  return(output)
}

zscore <- function(values) {
  return((values-mean(values))/sd(values))  
}

scalePalette <-function(color.palette, values, mid.value=0, resolution=100) {
  max.value = max(values)
  min.value = min(values)
  max.distance = abs(max.value-mid.value)
  min.distance = abs(min.value-mid.value)
  color.scale = color.palette(resolution)
  if (max.distance > min.distance) { 
    ratio = (min.distance+max.distance)/(2*max.distance)
    limit = ceiling(resolution*(1-ratio))
    color.scale = color.scale[limit:resolution]
  }
  else if (max.distance < min.distance) { 
    ratio = (min.distance+max.distance)/(2*min.distance)
    limit = floor(resolution*(ratio))
    color.scale = color.scale[0:limit] }
  return(color.scale)
}


print("Loaded colors...")

# Load example WW epistasis data:
test = FALSE
if (test) {
  data.colors = read.table("../example/mapepistasis_ww_fitness_ns_product_scores.mtx", header=TRUE)
  data.report = read.table("../example/mapepistasis_ww_fitness_ns_product_report.txt", header=TRUE)
  names(data.report)
  library(gplots)

  # Define color testing function:
  printColors <- function(colors=names(color.palettes), xlab="Epistasis", ylab="Fitness") {
    x = data.colors$epistasis
    y = log2(data.colors$fitness.fraction)
    matrix.colors = xtabs(positive.fraction - negative.fraction ~ i + j, data=data.report)
    for (color in colors) {
      print(paste("Processing:", color))
      dcols = densCols(x, y, colramp=color.palettes[[color]])
      plot(x, y, font.lab=2, col=dcols, pch=19, xlab=xlab, ylab=ylab, main=color)
      heatmap.2(matrix.colors, Rowv=FALSE, Colv=FALSE, col=color.palettes[[color]], scale="none", trace="none", main=color, symbreaks=TRUE, breaks=seq(-0.5, 0.5, by=0.05))
    }
  }

  pdf("graphs/mapCode-colors.pdf", width=5, height=5)
  printColors()
  dev.off()
  }


