## Dañando el archivo de Willi como un test.
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 80,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Results from 1 Provider:\n",
      "\n",
      "14 Results from the JSOCClient:\n",
      "         T_REC          TELESCOP  INSTRUME WAVELNTH CAR_ROT\n",
      "         str23            str7      str9   float64   int64 \n",
      "----------------------- -------- --------- -------- -------\n",
      "2013.11.06_00:00:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:02:15_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:04:30_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:06:45_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:09:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:11:15_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:13:30_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:15:45_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:18:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:20:15_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:22:30_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:24:45_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:27:00_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "2013.11.06_00:29:15_TAI  SDO/HMI HMI_SIDE1   6173.0    2143\n",
      "\n",
      "\n",
      "Export request pending. [id=\"JSOC_20200211_151_X_IN\", status=2]\n",
      "Waiting for 0 seconds...\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\r",
      "Files Downloaded:   0%|          | 0/14 [00:00<?, ?file/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "14 URLs found for download. Full request totalling 226MB\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n",
      "hmi.b_135s.20131106_000000_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "hmi.b_135s.20131106_000215_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000430_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000645_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000900_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "hmi.b_135s.20131106_001330_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n",
      "hmi.b_135s.20131106_001115_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_001800_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_001545_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_002015_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_002915_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "hmi.b_135s.20131106_002445_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_002700_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "Files Downloaded:   0%|          | 0/14 [00:01<?, ?file/s]      | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n"
     ]
    }
   ],
   "source": [
    "import astropy.units as u\n",
    "import matplotlib.pyplot as plt\n",
    "import drms\n",
    "import astropy.time\n",
    "c = drms.Client()  \n",
    "\n",
    "res = Fido.search(a.Time('2013-02-17T00:00:00', '2013-02-17T23:59:00'), a.jsoc.Series('hmi.b_135s'),\n",
    "                  a.jsoc.Segment('inclination'),\n",
    "                  a.jsoc.Notify('wjcaicedot@unal.edu.co'))\n",
    "res1 = Fido.search(a.Time('2013-11-06T00:00:00', '2013-11-06T23:59:00'), a.jsoc.Series('hmi.b_135s'),\n",
    "                  a.jsoc.Segment('inclination'),\n",
    "                  a.jsoc.Notify('wjcaicedot@unal.edu.co')) \n",
    "res2 = Fido.search(a.Time('2013-11-07T00:00:00', '2013-11-07T23:59:00'), a.jsoc.Series('hmi.b_135s'),\n",
    "                  a.jsoc.Notify('wjcaicedot@unal.edu.co')) \n",
    "res3 = Fido.search(a.Time('2013-11-08T00:00:00', '2013-11-08T23:59:00'), a.jsoc.Series('hmi.b_135s'),\n",
    "                  a.jsoc.Notify('wjcaicedot@unal.edu.co')) \n",
    "\n",
    "print(res1)\n",
    "downloaded_files = Fido.fetch(res1)  \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 84,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "<ExportRequest id=\"JSOC_20200211_151_X_IN\", status=2>\n",
      "JSOC_20200211_151_X_IN\n",
      "2\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\r",
      "Files Downloaded:   0%|          | 0/14 [00:00<?, ?file/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "14 URLs found for download. Full request totalling 226MB\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n",
      "hmi.b_135s.20131106_000000_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "hmi.b_135s.20131106_000215_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000430_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000900_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_000645_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "hmi.b_135s.20131106_001330_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_001545_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "hmi.b_135s.20131106_001115_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_002015_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_001800_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\u001b[A\n",
      "hmi.b_135s.20131106_002700_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\n",
      "\n",
      "\n",
      "hmi.b_135s.20131106_002445_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\n",
      "\n",
      "hmi.b_135s.20131106_002230_TAI.inclination.fits:   0%|          | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\n",
      "\n",
      "\n",
      "\n",
      "Files Downloaded:   0%|          | 0/14 [00:01<?, ?file/s]      | 0.00/16.9M [00:00<?, ?B/s]\u001b[A\u001b[A\u001b[A\u001b[A\n"
     ]
    }
   ],
   "source": [
    "from sunpy.net import jsoc\n",
    "client = jsoc.JSOCClient()\n",
    "\n",
    "\n",
    "res = client.search(a.Time('2013-11-06T00:00:00', '2013-11-06T23:59:00'), a.jsoc.Series('hmi.b_135s'),\n",
    "                  a.jsoc.Segment('inclination'),\n",
    "                  a.jsoc.Notify('wjcaicedot@unal.edu.co'))\n",
    "\n",
    "\n",
    "requests = client.request_data(res)  \n",
    "print(requests)\n",
    "print(requests.id)  \n",
    "print(requests.status)\n",
    "res = client.get_request(requests) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
