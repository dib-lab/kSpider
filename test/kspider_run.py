import _kSpider_internal as ks

ks.sourmash_sigs_indexing(sigs_dir = "/home/mabuelanin/dib-dev/kSpider_dev/test/sigs", kSize = 25)
ks.pairwise(index_prefix = "sigs", user_threads = 1)
# print(dir(ks))