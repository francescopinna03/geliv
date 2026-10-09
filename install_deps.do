foreach p in ivreg2 ranktest avar weakiv bcuse {
    ssc install `p', replace
}
