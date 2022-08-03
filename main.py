import math

FUEL_MODEL_G = [[2.5, 2, 5, 12, 0.5, 0.5], [2000, 109, 30, 8, 2000, 1500], 1, 25, 8000]
# index [0] - fuel loading list [1-hour, 10-hour, 100-hour, 1000-hour, herb, wood], (tons/acre)
# index [1] - SAVR list [1-hour, 10-hour, 100-hour, 1000-hour, herb, wood], (ft^-1)
# index [2] effective fuel bed depth, (ft)
# index [3] dead fuel moisture of extinction (%)
# index [4] dead and live fuel heat of combustion (Btu/lb)


def calculate_emc(temp, rh):
    """
    This function calculates equilibrium moisture content (EMC)
    :param temp: dry bulb temperature (°Fahrenheit)
    :param rh: relative humidity (%)
    :return: emc: equilibrium moisture content (%)
    """
    if rh < 10:
        emc = 0.03229 + 0.281073 * rh - 0.000578 * temp * rh
    elif 10 <= rh < 50:
        emc = 2.22749 + 0.160107 * rh - 0.014784 * temp
    else:
        emc = 21.0606 + 0.005565 * (rh ** 2) - 0.00035 * rh * temp - 0.483199 * rh
    return emc


# TODO
def calculate_potential_solar(lat, lon, elev, day):
    """
    This function calculates potential maximum solar insolation (maxsolar)
    :param lat: latitude (°)
    :param lon: longitude (°)
    :param elev: elevation (ft)
    :param day: !!!

    """
    sc = 1367  # solar constant (W/m^2)
    delta = 23.45 * math.sin(math.radians(360 / 365 * (284 + day)))  # declination angle (°)
    return delta


def calculate_sow(rad, ppt, maxrad):
    """
    This function determines state-of-weather code
    :param rad:
    :param ppt: average hourly precipitation in a day (in/hour)
    :param maxrad:
    :return: sow: state-of-weather code
    """
    percent_cloud = rad / maxrad
    if ppt < 0.01:
        if percent_cloud >= 0.91:
            sow = 0
        elif 0.73 <= percent_cloud < 0.91:
            sow = 1
        elif 0.50 <= percent_cloud < 0.73:
            sow = 2
        else:
            sow = 3
    elif 0.01 <= ppt < 0.05:
        sow = 5
    else:
        sow = 6
    return sow


def estimate_rh_temp_interface(temp, rh, sow):
    """
    This function estimates the relative humidity and dry bulb temperature of the air in immediate contact with the fuel
     elements
    :param temp: dry bulb temperature (°Fahrenheit)
    :param rh: relative humidity (%)
    :param sow: state-of-weather code
    :return: tempprm, rhprm
    """
    tempprm = None
    rhprm = None
    if sow == 0:
        tempprm = temp + 25
        rhprm = rh * 0.75
    elif sow == 1:
        tempprm = temp + 19
        rhprm = rh * 0.83
    elif sow == 2:
        tempprm = temp + 12
        rhprm = rh * 0.92
    elif sow == 3:
        tempprm = temp + 5
        rhprm = rh * 1
    return [tempprm, rhprm]


def calculate_mc1(emcprm, sow):
    """
    This functions calculates 1-hour timelag percent fuel moisture content (MC1)
    :param emcprm: boundary layer equilibrium moisture content (%)
    :param sow: state-of-weather code
    :return: mc1: calculated 1-hour timelag percent fuel moisture content
    """
    if sow > 4:
        mc1 = 35
    else:
        mc1 = 1.03 * emcprm
    return mc1


def calculate_mc10(emcprm, sow):
    """
    This functions calculates 10-hour timelag percent fuel moisture content (MC1)
    :param emcprm: boundary layer equilibrium moisture content (%)
    :param sow: state-of-weather code
    :return: mc10: calculated 10-hour timelag percent fuel moisture content
    """
    if sow > 4:
        mc10 = 35
    else:
        mc10 = 1.28 * emcprm
    return mc10


def calculate_mc10_model88(emcprm, mc100):
    """
    This functions calculates 10-hour timelag percent fuel moisture content (MC1)
    :param emcprm: boundary layer equilibrium moisture content (%)
    :param sow: state-of-weather code
    :return: mc10: calculated 10-hour timelag percent fuel moisture content
    """
    mc10 = 1.03 * emcprm * 0.8 + mc100 * 0.2
    return mc10


def calculate_daylight_hours(jdate, lat):
    """

    :param jdate:
    :param lat:
    :return:
    """
    converter = math.pi / 180
    decl = 0.41008 * math.sin((jdate - 82) * converter)  # solar declination (radians)
    phi = lat * converter  # station latitude (radians)
    daylit = 24 * (1 - math.acos(math.tan(phi) * math.tan(decl)) / math.pi)  # duration of daylight (hours)
    return daylit


def calculate_erc(mc1, mc10, mc100, mc1000, mcherb, mcwood, fuel_model):
    converter = 2000 / 43560 # too precise for us :) ---> but can be used for sure
    # converter = 0.046

    # Fuel loadings (lb/ft^2)
    w1 = fuel_model[0][0] * converter
    w10 = fuel_model[0][1] * converter
    w100 = fuel_model[0][2] * converter
    w1000 = fuel_model[0][3] * converter
    wherb = fuel_model[0][4] * converter
    wwood = fuel_model[0][5] * converter

    # SAVR of each fuel class (ft^-1)
    sg1 = fuel_model[1][0]
    sg10 = fuel_model[1][1]
    sg100 = fuel_model[1][2]
    sg1000 = fuel_model[1][3]
    sgherb = fuel_model[1][4]
    sgwood = fuel_model[1][5]

    fctcur = 4 / 3 - 1 / 90 * mcherb  # fraction of the herbaceous fuel loading transferred to the 1-hour class

    if fctcur > 1:
        fctcur = 1
    elif fctcur < 0:
        fctcur = 0

    wherbc = fctcur * wherb  # amount of the herbaceous fuel loading transferred to the 1-hour class
    w1p = w1 + wherbc  # 1-hour fuel loading and transferred herbaceous loading
    wherbp = wherb - wherbc  # amount of herbaceous loading left after transfer to 1-hour class

    # Fraction of the live and dead fuels made up of inert, noncombustible minerals
    std = 0.0555
    stl = 0.0555

    # Total live and dead fuel loading
    wtotd = w1p + w10 + w100 + w1000
    wtotl = wherbp + wwood
    wtot = wtotd + wtotl  # total fuel loading

    rhol = 32  # live fuel particle density (lb/ft^3)
    rhod = 32  # dead fuel particle density (lb/ft^3)
    depth = fuel_model[2]

    # Bulk density of the fuel bed
    rhobed = (wtot - w1000) / depth

    # Weighted fuel density
    rhobar = ((wtotl * rhol) + (wtotd * rhod)) / wtot

    # Packing ratio
    betbar = rhobed / rhobar

    # Weighting factors of each fuel class by fuel loading
    f1e = w1p / wtotd
    f10e = w10 / wtotd
    f100e = w100 / wtotd
    f1000e = w1000 / wtotd
    fherbe = wherbp / wtotl
    fwoode = wwood / wtotl

    # Weighting factors of dead and live fuels by fuel loading
    fdeade = wtotd / wtot
    flivee = wtotl / wtot

    # Net loadings of dead and live fuels by fuel loading
    wdedne = wtotd * (1 - std)
    wlivne = wtotl * (1 - stl)

    # Dead and live fuel characteristic surface-area-to-volume ratios by fuel loading
    sgbrde = (f1e * sg1) + (f10e * sg10) + (f100e * sg100) + (f1000e * sg1000)
    sgbrle = (fwoode * sgwood) + (fherbe * sgherb)

    # Characteristic surface-area-to-volume ratio by fuel loading
    sgbrte = (fdeade * sgbrde) + (flivee * sgbrle)

    # Optimum packing ratio by fuel loading
    betope = 3.348 * sgbrte ** (-0.8189)

    # Maximum reaction velocity by fuel loading
    gmamxe = (sgbrte ** 1.5) / (495 + 0.0594 * sgbrte ** 1.5)

    ade = 133 * sgbrte ** (-0.7913)
    # Optimum reaction velocity by fuel loading
    gmaope = gmamxe * (betbar / betope) ** ade * math.exp(ade * (1 - betbar / betope))

    # Weighted moisture contents of dead and live fuels by fuel loading
    wtmcde = (f1e * mc1) + (f10e * mc10) + (f100e * mc100) + (f1000e * mc1000)
    wtmcle = (fwoode * mcwood) + (fherbe * mcherb)

    # Moisture of extinction of dead fuels
    mxd = fuel_model[3]

    # Net fuel loading for each fuel class (lb/ft^2)
    w1n = w1p * (1 - std)
    w10n = w10 * (1 - std)
    w100n = w100 * (1 - std)
    wherbn = wherbp * (1 - stl)
    wwoodn = wwood * (1 - stl)

    # Heating numbers of ech fuel class
    hn1 = w1n * math.exp(-138 / sg1)
    hn10 = w10n * math.exp(-138 / sg10)
    hn100 = w100n * math.exp(-138 / sg100)
    hnherb = wherbn * math.exp(-500 / sgherb)
    hnwood = wwoodn * math.exp(-500 / sgwood)

    # Ratio of dead-to-live fuel heating numbers
    wrat = (hn1 + hn10 + hn100) / (hnherb + hnwood)

    # Weighted dead-fuel moisture content for live-fuel extinction moisture
    mclfe = ((mc1 * hn1) + (mc10 * hn10) + (mc100 * hn100)) / (hn1 + hn10 + hn100)

    # Moisture of extinction of live fuels
    mxl = (2.9 * wrat * (1 - mclfe / mxd) - 0.226) * 100
    if mxl < mxd:
        raise ValueError

    dedrte = wtmcde / mxd
    livrte = wtmcle / mxl

    # Moisture damping coefficients of dead and live fuels
    etamde = 1 - 2 * dedrte + 1.5 * dedrte ** 2 - 0.5 * dedrte ** 3
    if etamde > 1 or etamde < 0:
        raise ValueError
    etamle = 1 - 2 * livrte + 1.5 * livrte ** 2 - 0.5 * livrte ** 3
    if etamle > 1 or etamle < 0:
        raise ValueError

    # Fractions of the dead and live fuels made up of silica-free, noncombustible minerals
    sd = 0.01
    sl = 0.01

    # Mineral damping coefficient of live and dead fuels
    etasd = 0.174 * sd ** (-0.19)
    etasl = 0.174 * sl ** (-0.19)

    hd = fuel_model[4]
    hl = fuel_model[4]

    # Reaction intensity
    ire = gmaope * ((fdeade * wdedne * hd * etasd * etamde) + (flivee * wlivne * hl * etasl * etamle))

    # Surface area of each fuel class (ft^2)
    sa1 = (w1p / rhod) * sg1
    sa10 = (w10 / rhod) * sg10
    sa100 = (w100 / rhod) * sg100
    saherb = (wherbp / rhol) * sgherb
    sawood = (wwood / rhol) * sgwood

    # Total surface area of dead and live fuels
    sadead = sa1 + sa10 + sa100
    salive = saherb + sawood

    # Weighting factors of each fuel class by surface area
    f1 = sa1 / sadead
    f10 = sa10 / sadead
    f100 = sa100 / sadead
    fherb = saherb / salive
    fwood = sawood / salive

    # Weighting factors of dead and live fuels by surface area
    fdead = sadead / (sadead + salive)
    flive = salive / (sadead + salive)

    # Dead and live fuel characteristic surface-area-to-volume ratios by surface area
    sgbrd = (f1 * sg1) + (f10 * sg10) + (f100 * sg100)
    sgbrl = (fwood * sgwood) + (fherb * sgherb)

    # Characteristic surface-area-to-volume ratio by surface area
    sgbrt = (fdead * sgbrd) + (flive * sgbrl)

    # Residence time of the flaming front
    tau = 384 / sgbrt

    # Energy release component
    erc = 0.04 * ire * tau

    return erc


def calculate_emcbar(emcmin, emcmax, daylit):
    emcbar = (daylit * emcmin + (24 - daylit) * emcmax) / 24
    return emcbar


def calculate_bndryh(emcbar, pptdur):
    bndryh = ((24 - pptdur) * emcbar + pptdur * (0.5 * pptdur + 41)) / 24
    return bndryh


def calculate_mc100(ymc100, bndryh):
    mc100 = ymc100 + (bndryh - ymc100) * (1 - 0.87 * math.exp(-0.24))
    return mc100


def calculate_bndryt(emcbar, pptdur):
    bndryt = ((24 - pptdur) * emcbar + pptdur * (2.7 * pptdur + 76)) / 24
    return bndryt


def calculate_bdybar(set_bndryt):
    bdybar = sum(set_bndryt) / len(set_bndryt)
    return bdybar


def calculate_mc1000(pm1000, bdybar):
    mc1000 = pm1000 + (bdybar - pm1000) * (1 - 0.82 * math.exp(-0.168))
    return mc1000


print(round(calculate_erc(2.7, 3.2, 5.2, 6.9, 2.7, 60.0, FUEL_MODEL_G), 1))

# temp = 39
# rh = 84
# sow = 0
# ymc100 = 17.1
pm1000 = 20
pptdur = 15

daylit = calculate_daylight_hours(1, 39.74)
emc_min = calculate_emc(78, 21)
emc_max = calculate_emc(54, 31)
# print(emc_max)
emcbar = calculate_emcbar(emc_min, emc_max, daylit)
bndryh = calculate_bndryh(emcbar, pptdur)
# mc100 = calculate_mc100(ymc100, bndryh)

bndryt = calculate_bndryt(emcbar, pptdur)
bdybar = calculate_bdybar([20, 20, 20, 20, 20, 20, bndryt])
mc1000 = calculate_mc1000(pm1000, bdybar)
# print(bndryt, mc1000)

# values = estimate_rh_temp_interface(temp, rh, sow)
# emc = calculate_emc(values[0], values[1])
# mc1 = calculate_mc1(emc, sow)
# print(mc1)
# mc10 = calculate_mc10_model88(emc, mc100)
# print(mc1, mc10, mc100)
