licenses(["notice"])

cc_library(
    name = "Dendrogram",
    hdrs = ["Dendrogram.h"],
    deps = [
        "@gbbs//gbbs",
        "//Dendrogram/common:leftist_heap",
        "//Dendrogram/common:skew_heap",
        "//Dendrogram/common:pairing_heap",
        "//Dendrogram/common:union_find",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [
        ":Dendrogram",
        "//utils:benchmark",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
