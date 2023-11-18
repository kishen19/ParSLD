#licenses(["notice"])
#
#cc_library(
#    name = "TreeContraction",
#    hdrs = ["TreeContraction.h"],
#    deps = [
#        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
#    ],
#)
#
#cc_library(
#    name = "Dendrogram",
#    hdrs = ["Dendrogram.h"],
#    deps = [
#        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
#        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
#        "@gbbs//gbbs",
#        "//Dendrogram/heaps:leftist_heap",
#        "//Dendrogram/heaps:skew_heap",
#        "//Dendrogram/heaps:pairing_heap",
#        ":TreeContraction",
#    ],
#)
#
#cc_binary(
#    name = "Dendrogram_main",
#    srcs = ["Dendrogram.cc"],
#    deps = [
#        ":Dendrogram",
#    ],
#)
#
#package(
#    default_visibility = ["//visibility:public"],
#)
