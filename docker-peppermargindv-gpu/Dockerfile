FROM kishwars/pepper_deepvariant:r0.8-gpu

# File Author / Maintainer
MAINTAINER David Spencer <dspencerlab@gmail.com>

LABEL Image with PEPPER-Margin-DeepVariant for ONT processing

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

ENV PATH=$PATH:/opt/margin_dir/build/:/opt/deepvariant/bin/
