
#Method 2
raagaCount = 0
for r in raagas:
    recCount = 0
    for mbid in raagaMBIDs[r][:5]:
        data = pickle.load(file(methodPath+"/"+mbid+".pickle"))
        [n, be] = np.histogram(data["G3"], bins=1200, density=True)
        bc = (be[:-1]+be[1:])/2.0
        ns = gaussian_filter(n, 7)
        ax[raagaCount].plot(bc, ns, color=colors[recCount], ls=styles[recCount], lw="1.5")
        recCount += 1
    raagaCount += 1
        

#Varnams
#raagas = ['begada', 'abhogi', 'mohanam', 'shree']
raagas = ["shree"]
rindex = 0
for x in [0, 1]:
    for y in [0, 1]:
        r = raagas[rindex]
        rindex += 1
        artists = listdir("/home/gopal/workspace/intonationLib/data/varnam-analysis/recorded/distributions/"+r)
        artists = [i[:-7] for i in artists]
        if "dharini" in artists: artists.remove("dharini") 
        for s in ["R"]:#distributions[r+"_"+artists[0]].keys():
            count = 0
            for a in artists:
                [n, binEdges] = histogram(distributions[r+"_"+a][s], bins=1200, density=True)
                nS = gaussian_filter(n, 5)
                binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
                ax.plot(binCenters, nS, color=colors[count], ls=styles[count], lw="1.5")
                #plot(binCenters, nS)
                count += 1
