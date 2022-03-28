# Multi sensor waterextraction
Advances in remote sensing technologies helped us to improve the monitoring of wetlands; however, detecting the presence of water under vegetation is still a challenge for correctly delineating the water extent. To address this issue and better detect the presence of water below vegetation, we employ different polarization of SAR sentinel-1 data in combination with optical sentinel-2. After preprocessing the images, we use the K-means clustering algorithm provided in the cloud computing platform of Google Earth Engine, to detect the increased backscatter coming from flooded vegetation duo to the double-bounce of the radar signal. We also take advantage of the high-resolution national land cover of Sweden as an ancillary layer to extract only the relevant information in our study area. Finally, we compare our results with hydroclimatic and field data gathered from the study area. Our workflow improves water-extent delineation in Swedish wetlands by 20% on average by detecting hidden water below the vegetation, which is generally not recognized by optical methods. The proposed method can be extended to monitor and study wetlands’ water availability and changes, contributing to the increase of their resilience to anthropogenic pressures and climate change.
