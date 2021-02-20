# ![goGPS](https://github.com/goGPS-Project/goGPS_graphics/blob/master/logo/PNGs/goGPS_logo_128.png?raw=true) goGPS 1.0 Open Edition

_Welcome to the goGPS wiki pages, here we will try to write useful information for users and developers._

**Disclaimer**: Some links might be broken in the wiki here on GitHub. If you need to access documentation or search for more information on the software visit our website at [gogps-project.github.io](https://gogps-project.github.io). 

**goGPS** is a software created for processing GNSS raw data. It was originally written specifically to work with GPS single-frequency low-cost receivers but now it can fully exploit **multi-constellation**, **multi-frequency**, **multi-tracking** observations. goGPS implements multiple algorithms to analyze the data, and at the moment these include two main Least Squares (LS) engines: one working on combination of observables (e.g. iono-free observations) and one able to use all the frequencies and trackings logged without performing any combinations (ionospheric delay are parameters of the normal equations).
The **combined and uncombined engines** can both compute Precise Point Positioning (PPP) solutions and Network adjustments (NET).

**Note:** At the moment the software is focused on the processing of permanent stations (geodetic or low-cost), it does not yet include the possibility to analyze moving receivers.

## Index of contents:

1. [Introduction](https://gogps-project.github.io/wiki/Introduction)
   - [A bit of history](https://gogps-project.github.io/wiki/Introduction#a-bit-of-history)
   - [About](https://gogps-project.github.io/wiki/Introduction#about)
2. [Installation](https://gogps-project.github.io/wiki/Installation)
   - [Requirements](https://gogps-project.github.io/wiki/Installation#requirements)
   - [Install procedure](https://gogps-project.github.io/wiki/Installation#installation)
   - [goGPS directories](https://gogps-project.github.io/wiki/Installation#directories-of-gogps)
   - [Execution](https://gogps-project.github.io/wiki/Installation#execution) 
3. [Settings](https://gogps-project.github.io/wiki/Settings)
   - [Index of contents](https://gogps-project.github.io/wiki/Settings#index-of-contents)
   - [Menu](https://gogps-project.github.io/wiki/Settings#menu)
   - [Sidebar](https://gogps-project.github.io/wiki/Settings#sidebar)
   - [Tabs](https://gogps-project.github.io/wiki/Settings#tabs)
      * [Advanced](https://gogps-project.github.io/wiki/Settings#tab-advanced)
      * [Resources](https://gogps-project.github.io/wiki/Settings#tab-resources)
      * [Commands](https://gogps-project.github.io/wiki/Settings#tab-commands)
      * [Data Sources](https://gogps-project.github.io/wiki/Settings#tab-data-sources)
      * [Receiver Info](https://gogps-project.github.io/wiki/Settings#tab-receiver-info)
      * [Pre-Processing](https://gogps-project.github.io/wiki/Settings#tab-pre-processing)
      * [Processing](https://gogps-project.github.io/wiki/Settings#tab-processing)
      * [Parametrization PPP](https://gogps-project.github.io/wiki/Settings#tab-parametrization-ppp)
      * [Parametrization Network](https://gogps-project.github.io/wiki/Settings#tab-parametrization-network)
      * [Atmosphere](https://gogps-project.github.io/wiki/Settings#tab-atmosphere)
   - [Bottom Bar](https://gogps-project.github.io/wiki/Settings#tab-bottom-bar)
4. [Command language](https://gogps-project.github.io/wiki/Command-language)
   * [Index of contents](https://gogps-project.github.io/wiki/Command-language#index-of-contents)
   * [Receiver data organization](https://gogps-project.github.io/wiki/Command-language#receiver-data-organization)
   * [Minimal script](https://gogps-project.github.io/wiki/Command-language#minimal-script)
      * [PUSHOUT](https://gogps-project.github.io/wiki/Command-language#pushout)
   * [How to select sessions](https://gogps-project.github.io/wiki/Command-language#how-to-select-sessions)
   * [How to select receivers](https://gogps-project.github.io/wiki/Command-language#how-to-select-receivers)
   * [Main commands](https://gogps-project.github.io/wiki/Command-language#main-commands)
      * [LOAD](https://gogps-project.github.io/wiki/Command-language#load)
      * [PREPRO](https://gogps-project.github.io/wiki/Command-language#prepro)
      * [PPP](https://gogps-project.github.io/wiki/Command-language#ppp)
      * [NET](https://gogps-project.github.io/wiki/Command-language#net)
   * [Looping on receivers](https://gogps-project.github.io/wiki/Command-language#looping-on-receivers)
   * [Parallel execution](https://gogps-project.github.io/wiki/Command-language#parallel-execution)
      * [PINIT](https://gogps-project.github.io/wiki/Command-language#pinit)
      * [PKILL](https://gogps-project.github.io/wiki/Command-language#pkill)
      * [PAR](https://gogps-project.github.io/wiki/Command-language#par)
   * [Plots and export](https://gogps-project.github.io/wiki/Command-language#plots-and-export)
      * [SHOW](https://gogps-project.github.io/wiki/Command-language#show)
      * [EXPORT](https://gogps-project.github.io/wiki/Command-language#export)
   * [Validation](https://gogps-project.github.io/wiki/Command-language#validation)
      * [VALIDATE](https://gogps-project.github.io/wiki/Command-language#validate)
   * [Ionospheric interpolation](https://gogps-project.github.io/wiki/Command-language#ionospheric-interpolation)
      * [SEID](https://gogps-project.github.io/wiki/Command-language#seid)
      * [SID](https://gogps-project.github.io/wiki/Command-language#sid)
      * [REMIONO](https://gogps-project.github.io/wiki/Command-language#remiono)
   * [Multipath management](https://gogps-project.github.io/wiki/Command-language#multipath-management)
      * [MPEST](https://gogps-project.github.io/wiki/Command-language#mpest)
   * [Other Commands](https://gogps-project.github.io/wiki/Command-language#other-commands)
      * [RENAME](https://gogps-project.github.io/wiki/Command-language#rename)
      * [EMPTY](https://gogps-project.github.io/wiki/Command-language#empty)
      * [EMPTYWORK](https://gogps-project.github.io/wiki/Command-language#emptywork)
      * [EMPTYOUT](https://gogps-project.github.io/wiki/Command-language#emptyout)
      * [AZEL](https://gogps-project.github.io/wiki/Command-language#azel)
      * [BASICPP](https://gogps-project.github.io/wiki/Command-language#basicpp)
      * [CODEPP](https://gogps-project.github.io/wiki/Command-language#codepp)
      * [OUTDET](https://gogps-project.github.io/wiki/Command-language#outdet)
      * [FIXPOS](https://gogps-project.github.io/wiki/Command-language#fixpos)
      * [KEEP](https://gogps-project.github.io/wiki/Command-language#keep)
      * [SYNC](https://gogps-project.github.io/wiki/Command-language#sync)
      * [REMSAT](https://gogps-project.github.io/wiki/Command-language#remsat)
      * [REMOBS](https://gogps-project.github.io/wiki/Command-language#remobs)
      * [REMTMP](https://gogps-project.github.io/wiki/Command-language#remtmp)
1. **[How to create a new project](https://gogps-project.github.io/wiki/How-to-create-a-new-project)**

### Coding goGPS
For a bit old presentation on the goGPS coding have a look at this presentation:
[Coding goGPS](https://github.com/goGPS-Project/goGPS_MATLAB/blob/goGPS_1.0_beta/docs/coding%20goGPS.pdf)
