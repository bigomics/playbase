test_that("example files exist", {
  all_files <- playbase::example_file()
  expect(
      sum(!c(
          'counts.csv',
          'contrasts.csv',
          'samples.csv'
      ) %in% all_files) == 0,
      'Cant find some external data files in playbase_example()'
  )
})
