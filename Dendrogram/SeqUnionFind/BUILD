licenses(["notice"])

cc_library(
    name = "Dendrogram",
    hdrs = ["Dendrogram.h"],
    deps = [
        "//Dendrogram/common:union_find",
        "@gbbs//gbbs",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [":Dendrogram"],
)

package(
    default_visibility = ["//visibility:public"],
)
