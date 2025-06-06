## ----setup, include=FALSE-----------------------------------
#AnnArbor 
#linkcolor: #191970
knitr::opts_chunk$set(echo = TRUE, comment = NA,cache=TRUE)
library(tidyverse)


## ---- echo=FALSE, out.width='10%',fig.align='right'---------
knitr::include_graphics("Imgs/hex-tidyverse.png")


## ---- echo=F, out.width='5%', fig.align='right'-------------
knitr::include_graphics("Imgs/hex-tidyverse.png")


## ----warning=FALSE, message=TRUE----------------------------
#install.packages("tidyverse")
library(tidyverse)


## ---- echo=F, out.width='60%', fig.align="center"-----------
knitr::include_graphics("Imgs/data-science.png")


## ---- echo=F, out.width='30%', fig.align="center"-----------
knitr::include_graphics("Imgs/data-science.png")


## ----serpiente_camello, echo=FALSE,fig.align='center',out.width="40%"----
knitr::include_graphics("Imgs/serpiente_camello.png")


## ----tidydata, echo=FALSE,fig.align='center',out.width="60%"----
knitr::include_graphics("Imgs/tidy_data.PNG")


## ---- warning=FALSE-----------------------------------------
#install.packages("palmerpenguins",dep=TRUE)
library("palmerpenguins")
print(penguins, width = 50)


## ---- echo=F, warning=F-------------------------------------
set.seed(123)

penguins %>% 
  group_by(species, island) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pivot_wider(names_from = island, values_from = n) %>% 
  unnest(cols = c(Biscoe, Dream, Torgersen))


## ---- echo=F------------------------------------------------
penguins %>% 
  select(species, island, sex, year) %>% 
  unite(col, species, sex) %>% 
  sample_n(5)


## ---- echo=F, message=F, warning=F--------------------------
penguins %>% 
  select(bill_length_mm, bill_depth_mm, flipper_length_mm) %>% 
  corrr::correlate(method = "pearson")


## ---- echo=F,message=F, warning=F---------------------------
penguins %>% 
  select(species, island, sex) %>% 
  sample_n(3) %>% 
  bind_rows(
    mtcars %>%
      tibble::rownames_to_column("model") %>% 
      select(model, mpg, cyl) %>% 
      sample_n(3)
  )


## ---- echo=F, message=F, warning=F--------------------------
penguins %>% 
  select(species, island, sex) %>% 
  sample_n(3) %>% 
  bind_rows(
    mtcars %>%
      tibble::rownames_to_column("model") %>% 
      select(model, mpg, cyl) %>% 
      sample_n(3)
  )


## ----plot_tidy, echo=F,fig.align="center",include=FALSE, out.width="30%"----
png("Imgs/plot_tidy.png")
set.seed(33)
x=rnorm(1000)
y=rnorm(1000)
colores=cut(x,4)
levels(colores)=c("red","blue","green","brown")
par(mfrow=c(1,2))
plot(x=x,y=y,pch=15,col=sample(c("red","blue","green","brown"),500,replace=TRUE),
axes=FALSE,xlab="",ylab="",main="Datos NO tidy")
plot(x=x,y=y,pch=15,col=as.character(colores),axes=FALSE,xlab="",ylab="",main="Datos tidy")
par(mfrow=c(1,1))
dev.off()


## ----tidy_data, echo=FALSE,fig.align='center',out.width="50%"----
knitr::include_graphics("Imgs/plot_tidy.png")


## ---- echo=F, fig.align='center', out.width='30%'-----------
knitr::include_graphics("Imgs/pinguinos_madagascar.jpg")


## ---- echo=F, fig.align='center', out.width='35%'-----------
knitr::include_graphics(c("Imgs/lter_penguins.png","Imgs/culmen_depth.png"))


## ----fig.align='right',out.width="5%",fig.align='right',echo=FALSE----
knitr::include_graphics("Imgs/logo_pipe.png")


## ---- eval=FALSE--------------------------------------------
## rnorm(200) %>%
## matrix(ncol = 2) %T>%
## plot %>% # plot no suele retornar nada
## colSums


## ---- message=FALSE, eval=FALSE-----------------------------
## library(magrittr)
## iris %>%
##   subset(Sepal.Length > mean(Sepal.Length)) %$%
##   cor(Sepal.Length, Sepal.Width)


## ----out.width="50%", fig.align='center'--------------------
data.frame(z = rnorm(100)) %$%  ts.plot(z)


## ---- eval=F------------------------------------------------
## x %>% f # equivalente a: f(x)
## x %>% f(y) # equivalente a: f(x, y)
## x %>% f %>% g %>% h # equivalente a: h(g(f(x)))


## ---- eval=F------------------------------------------------
## x %>% f(.) # equivalente a: x %>% f
## x %>% f(y, .) # equivalente a: f(y, x)
## x %>% f(y, z = .) # equivalente a: f(y, z = x)
## x %>% f(y = nrow(.),
##         z = ncol(.))  # equivalente a: f(x, y = nrow(x), z = ncol(x))


## ---- eval=F------------------------------------------------
## f <- . %>% cos %>% sin # equivalente a: f <- function(.) sin(cos(.))

## ---- eval=F------------------------------------------------
## f(20) # equivalente a: la tubería 20 %>% cos %>% sin


## ---- eval=F------------------------------------------------
## mean(subset(penguins, year == 2007)$body_mass_g, na.rm = T)
## 
## # alternativamente:
## peng_bm_2007 <- subset(penguins, año == 2007)$body_mass_g
## media(peng_bm_2007, na.rm = T)


## ---- eval=F------------------------------------------------
## penguins %>%
##   subset(año == 2007) %>%
##   .$body_mass_g %>%
##   mean(na.rm = T)


## -----------------------------------------------------------
mtcars |> head()  #  es lo mismo que head(mtcars)
mtcars |> head(2) #  es lo mismo que  head(mtcars, 2)
mtcars |> subset(cyl == 4) |> nrow()  


## -----------------------------------------------------------
penguins[1:5, c("island", "bill_length_mm" )] %T>% 
  print %>% .$"bill_length_mm"  %>%
  mean(na.rm=T)


## ---- eval=FALSE--------------------------------------------
## plot(penguins$species, penguins$bill_length_mm)


## ---- fig.width=6, fig.asp=0.618, fig.retina=3,out.width="50%", fig.align='center'----
penguins %$% 
  plot(species,bill_length_mm) 
# 


## ---- results='asis'----------------------------------------
variable <- penguins$bill_length_mm
variable %<>% mean(na.rm=T)
variable


## -----------------------------------------------------------
tibble(
  x = c("a", "b"),
  y = c(1, 2),
  z = c(T, F)
)


## -----------------------------------------------------------
tribble(
  ~x, ~y, ~z,
  "a", 1, T,
  "b", 2, F
)


## -----------------------------------------------------------
data.frame(
  x = c("a", "b"),
  y = c(1, 2),
  z = c(T, F)
) %>% 
as_tibble


## -----------------------------------------------------------
c(x = "a", y = "b", z = 1) %>%
  enframe(name = "x", value = "y")


## -----------------------------------------------------------
penguins


## -----------------------------------------------------------
data.frame(penguins)


## -----------------------------------------------------------
penguins %>% glimpse


## -----------------------------------------------------------
data.frame(penguins) %>% .[, "species"] %>% class


## -----------------------------------------------------------
penguins[, "species"] %>% class


## -----------------------------------------------------------
names(data.frame(penguins))
head(data.frame(penguins)$spec)


## -----------------------------------------------------------
head(penguins$spec)


## -----------------------------------------------------------
head(penguins$species)


## -----------------------------------------------------------
data <- read_csv(file = "data/penguins.csv")


## -----------------------------------------------------------
data <- read_csv(file = "data/penguins.csv", col_select = c(species, island))


## -----------------------------------------------------------
data <- read_csv(file = "./data/penguins.csv",
                 col_names = paste("Var", 1:8, sep = "_"))


## -----------------------------------------------------------
data <- read_csv(file = "./data/penguins.csv", skip = 5)


## ---- eval=F------------------------------------------------
## read_csv(
##   archivo = "./data/penguins.csv",
##   col_types = cols(
##     species = col_character(),
##     año = col_datetime(formato = "%Y"),
##     isla = col_skip())
##   )


## ---- eval=F------------------------------------------------
## read_csv(file = "./data/penguins.csv", guess_max = 2000)


## -----------------------------------------------------------
penguins %>% 
  write_rds(file = "./data/penguins.rds")


## -----------------------------------------------------------
penguins <- read_rds(file = "./data/penguins.rds")


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
penguins


## ---- out.width="50%",fig.align='center',echo=FALSE---------
knitr::include_graphics("Imgs/pivotting.png") 


## -----------------------------------------------------------
long_penguins <- penguins %>% 
  pivot_longer(
    cols = c(species, island),
    names_to = "variable", values_to = "valor"
  )

long_penguins %>% glimpse



## -----------------------------------------------------------
long_penguins %>% 
  pivot_wider(
    names_from = "variable", values_from = "valor"
  ) %>%
glimpse


## -----------------------------------------------------------
nested_penguins <- penguins %>% 
    nest(nested_data = 
           c(island, bill_length_mm, 
             bill_depth_mm,flipper_length_mm,
             body_mass_g, sex))
nested_penguins


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
nested_penguins %>% purrr::pluck("nested_data", 1)


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
nested_penguins %>% unnest(cols = c(nested_data)) 


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png")


## -----------------------------------------------------------
nested_penguins %>% hoist(nested_data, hoisted_col = "bill_length_mm")


## -----------------------------------------------------------
penguins %>% unite(col = "specie_sex",
                   c(species, sex), sep = "_", remove = T)


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
penguins %>% separate(bill_length_mm, sep = 2, into = c("cm", "mm"))


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
penguins %>% separate_rows(island, sep = "s", convert = T)


## ---- echo=F------------------------------------------------
incompl_penguins <- tibble(
  species = c(rep("Adelie", 2), rep("Gentoo", 1), rep("Chinstrap", 1)),
  year = c(2007, 2008, 2008, 2007),
  measurement = c(rnorm(3, mean = 50, sd = 15), NA)
)


## -----------------------------------------------------------
incompl_penguins


## -----------------------------------------------------------
incompl_penguins %>% 
  complete(species, year, fill = list(measurement = NA))


## -----------------------------------------------------------
incompl_penguins %>% 
  drop_na(measurement)


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
incompl_penguins %>% 
  fill(measurement, .direction = "down")


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
incompl_penguins %>%
  replace_na(replace = list(measurement = mean(.$measurement, na.rm = T)))


## ---- out.width="5%",fig.align='right',echo=FALSE-----------
knitr::include_graphics("Imgs/logo_tidyr.png") 


## -----------------------------------------------------------
penguins %>% 
  filter(species == "Adelie")


## -----------------------------------------------------------
penguins %>% 
  filter(is.na(bill_length_mm) == T)
  # quitar los Nas
  #filter(!is.na(bill_length_mm) == F)


## -----------------------------------------------------------
penguins %>% 
  filter(between(body_mass_g, 3800, 4000) & (year < 2008 | year > 2008))


## -----------------------------------------------------------
penguins %>% 
  slice(23:27)


## -----------------------------------------------------------
penguins %>% 
  slice_head(n = 5) 
# alternativamente: 
#slice_head(frac = 0.05)


## -----------------------------------------------------------
penguins %>% 
  slice_sample(n = 5)


## -----------------------------------------------------------
penguins %>% 
  slice_max(bill_length_mm, n = 5)


## -----------------------------------------------------------
penguins %>% 
  arrange(body_mass_g) %>% 
  slice_head(n = 5)  # equivalentente a: slice_min(body_mass_g, n = 3)


## ---- out.width="10%",fig.align='right',echo=FALSE----------
knitr::include_graphics("Imgs/logo_dplyr.PNG") 


## -----------------------------------------------------------
penguins %>% 
  select(1:3) %>% 
  glimpse


## -----------------------------------------------------------
penguins %>% 
  select(species, island, bill_length_mm) %>% 
  glimpse


## -----------------------------------------------------------
penguins %>% 
  select(everything()) %>% 
  glimpse

# select(last_col())


## -----------------------------------------------------------
penguins %>% 
  select(starts_with("bill")) %>% 
  glimpse
# ends_with()


## -----------------------------------------------------------
penguins %>% 
  select(contains("e") & contains("a")) %>% 
  glimpse


## -----------------------------------------------------------
penguins %>% 
  select(matches("_\\w*_mm$")) %>% 
  glimpse


## -----------------------------------------------------------
penguins %>% 
  select(where(is.numeric)) %>% 
  glimpse


## ---- eval=F------------------------------------------------
## penguins %>%
##   select(ends_with("mm"))


## ---- eval=F------------------------------------------------
## penguins %>%
##   select(-contains("mm"))


## ---- eval=F------------------------------------------------
## penguins %>%
##   select(where(~ is.numeric(.))) %>%  # select(where(is.numeric))
##   select(where(~ mean(., na.rm = T) > 1000))


## -----------------------------------------------------------
penguins %>% rename(bm = body_mass_g, gender = sex) %>% 
  colnames()


## -----------------------------------------------------------
penguins %>% rename_with(.fn = toupper, .cols = contains("mm")) %>% 
  colnames()


## -----------------------------------------------------------
penguins %>% 
  relocate(species, .after = body_mass_g) %>%
  relocate(sex, .before = species) %>%
  relocate(island, .after = last_col()) %>%
  colnames()


## -----------------------------------------------------------
penguins %>% 
  mutate(bm_kg = body_mass_g / 1000, .keep = "all", .after = island) %>% 
  slice_head(n = 5)


## -----------------------------------------------------------
penguins %>% 
  mutate(
    sex_binary = case_when(
      sex == "male" ~ 1,
      sex == "female" ~ 0),
    .keep = "all", .after = island
  ) %>% 
  slice_head(n = 3)


## -----------------------------------------------------------
penguins %>% 
  mutate(
    across(contains("mm"), ~ . / 1000),
    .keep = "all"
  ) %>% 
  slice_head(n = 3)


## -----------------------------------------------------------
penguins %>% 
  mutate(
    across(where(is.character), as.factor),
    .keep = "all"
  ) %>% 
  slice_head(n = 3)


## ---- out.width="10%",fig.align='right',echo=FALSE----------
knitr::include_graphics("Imgs/logo_dplyr.PNG") 


## -----------------------------------------------------------
penguins %>% group_by(species)


## -----------------------------------------------------------
penguins %>% group_by(species) %>% 
  summarise(count = n(), .groups = "drop")


## -----------------------------------------------------------
penguins %>% group_by(species, sex) %>% summarise(count = n(), .groups = "drop")


## -----------------------------------------------------------
penguins %>%
  group_by(species) %>%
  summarise(
    across(contains("mm"), ~ mean(., na.rm = T), .names = "{.col}_avg"),
    .groups = "drop"
  )


## -----------------------------------------------------------
penguins %>% 
  group_by(species) %>% 
  group_by(year, .add = T)   # equivalente a: group_by(species, year)


## -----------------------------------------------------------
penguins %>%
  group_by(species) %>%
  summarise(
    across(
      contains("mm"),
      list(avg = ~ mean(., na.rm = T), sd = ~ sd(., na.rm = T)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


## -----------------------------------------------------------
penguins %>%
  group_by(species) %>% 
  mutate(stand_bm = (body_mass_g - mean(body_mass_g, na.rm = T))
         / sd(body_mass_g, na.rm = T)) %>% 
  glimpse


## ---- eval=FALSE--------------------------------------------
## bm_breaks <- mean(penguins$body_mass_g,
##                   na.rm = T) - (-3:3) *
##   sd(penguins$body_mass_g,na.rm = T)
## 
## penguins %>%
##   group_by(species, bm_bin = cut(body_mass_g, breaks = bm_breaks)) %>%
##   summarise(count = n(), .groups = "drop")


## ---- echo=FALSE--------------------------------------------
bm_breaks <- mean(penguins$body_mass_g, 
                  na.rm = T) - (-3:3) *
  sd(penguins$body_mass_g,na.rm = T)

penguins %>% 
  group_by(species, bm_bin = cut(body_mass_g, breaks = bm_breaks)) %>%
  summarise(count = n(), .groups = "drop")


## -----------------------------------------------------------
penguins %>% 
  group_by(species, island) %>% 
  filter(flipper_length_mm == max(flipper_length_mm, na.rm = T))


## -----------------------------------------------------------
penguins %>% 
  group_by(species, year) %>% 
  tidyr::nest()


## -----------------------------------------------------------
penguins %>% 
  distinct(species, island)


## ---- eval=FALSE--------------------------------------------
## penguins %>%
##   pull(year)  # equivalente a: penguins$year


## -----------------------------------------------------------
penguins %>% select(species, island, body_mass_g) %>% 
  mutate(penguin_size = if_else(body_mass_g < 3500,
                                "tiny penguin",
                                "big penguin"))


## -----------------------------------------------------------
penguins %>% select(species, body_mass_g) %>% 
  mutate(lagged_bm = lag(body_mass_g, n = 1))


## ---- out.width="40%",fig.align='center',echo=FALSE---------
knitr::include_graphics("Imgs/combinar_tablas.PNG") 


## ---- echo=F, out.height='30%', out.width='30%', fig.align='center'----
knitr::include_graphics("Imgs/ggplot2.PNG")


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE----
penguins %>% 
  ggplot(aes(x=bill_length_mm, y = flipper_length_mm)) +
  geom_point(na.rm = TRUE)


## ---- echo=F, out.height='60%', out.width='70%', fig.align='center'----
knitr::include_graphics("Imgs/escalas.PNG")


## ---- echo=F, out.height='60%', out.width='70%', fig.align='center'----
knitr::include_graphics("Imgs/stats.PNG")


## ---- echo=F, out.height='50%', out.width='40%', fig.align='center'----
knitr::include_graphics("Imgs/coordenadas.PNG")


## ---- echo=F, out.height='50%', out.width='40%', fig.align='center'----
knitr::include_graphics("Imgs/facetas.PNG")


## ---- echo=F, out.height='50%', out.width='40%', fig.align='center'----
knitr::include_graphics("Imgs/posicion.PNG")


## ---- echo=F, out.height='40%', out.width='40%', fig.align='center'----
knitr::include_graphics("Imgs/temas.PNG")


## ---- eval=FALSE--------------------------------------------
## penguins %>%
##   ggplot(aes(x = species)) +
##   geom_bar(fill="blue") +
##   labs(x="Especie", y="Número de pingüinos") +
##   theme_bw() +
##   theme(axis.text = element_text(size=20),
##         axis.title = element_text(size=20, face = "bold"))


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, echo=FALSE----
penguins %>% 
  ggplot(aes(x = species)) +
  geom_bar(fill="blue") + 
  labs(x="Especie", y="Número de pingüinos") +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"))
  


## ---- eval=FALSE--------------------------------------------
## penguins %>% ggplot() +
##   geom_bar(aes(species, fill=island),
##            position="dodge") + coord_flip() +
##   guides(fill = guide_legend(title = "Isla")) +
##   labs(x="Número de pingüinos", y="Especie") +
##   theme_bw() +
##   theme(axis.text = element_text(size=20),
##         axis.title = element_text(size=20, face = "bold"),
##         legend.title = element_text(size=20))
## 


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, echo=FALSE----
penguins %>% ggplot() + 
  geom_bar(aes(species, fill=island),
           position="dodge") + coord_flip() +
  guides(fill = guide_legend(title = "Isla")) +
  labs(y="Número de pingüinos", x="Especie") +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20)) 
  


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, echo=FALSE----
penguins %>% ggplot() + 
  geom_bar(aes(species, fill=island),
           position="fill") + coord_flip() +
  guides(fill = guide_legend(title = "Isla")) +
  labs(y="Proporción de pingüinos", x="Especie") +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20))


## ---- out.height='40%', out.width='60%', fig.align='center', message=FALSE----
penguins %>% 
  ggplot(aes(x = flipper_length_mm)) +
  geom_histogram(na.rm = TRUE) +
    labs(x="Longitud de la aleta en mm", 
         y="Frecuencia absoluta") + 
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"))


## ---- eval=FALSE--------------------------------------------
## ggplot(data = penguins, aes(x = flipper_length_mm)) +
##   geom_histogram(aes(fill = species),
##                  alpha = 0.5,
##                  position = "identity",
##                  na.rm = TRUE) +
##   scale_fill_manual(values = c("darkorange","purple","cyan4")) +
##   labs(x = "Longitud de la aleta en mm",
##        y = "Frecuencia absoluta") +
##   guides(fill = guide_legend(title = "Especie")) +
##   theme_bw() +
##   theme(axis.text = element_text(size=20),
##         axis.title = element_text(size=20, face = "bold"),
##         legend.title = element_text(size=20))


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, echo=FALSE----
ggplot(data = penguins, aes(x = flipper_length_mm)) +
  geom_histogram(aes(fill = species), 
                 alpha = 0.5, 
                 position = "identity",
                 na.rm = TRUE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "Longitud de la aleta en mm",
       y = "Frecuencia absoluta") +
  guides(fill = guide_legend(title = "Especie")) +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20))


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, echo=FALSE, message=FALSE----
ggplot(data = penguins, aes(x = flipper_length_mm)) +
  geom_histogram(aes(fill = species), 
                 alpha = 0.5, 
                 position = "identity",
                 na.rm = TRUE) +
  facet_grid(.~species) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "Longitud de la aleta en mm",
       y = "Frecuencia absoluta") +
  guides(fill = guide_legend(title = "Especie")) +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20))


## ---- eval=FALSE--------------------------------------------
## ggplot(data = penguins, aes(x = species, y = flipper_length_mm)) +
##   geom_boxplot(aes(color = species), width = 0.3,
##                show.legend = FALSE) +
##   geom_jitter(aes(color = species), alpha = 0.5,
##               show.legend = FALSE,
##               position = position_jitter(width = 0.2, seed = 0)) +
##   scale_color_manual(values = c("darkorange","purple","cyan4")) +
##   labs(x = "Epecie", y = "Longitud de la aleta en mm")


## ---- out.height='50%', out.width='60%', fig.align='center', message=FALSE, warning=FALSE, echo=FALSE----
ggplot(data = penguins, aes(x = species, 
                            y = flipper_length_mm)) +
  geom_boxplot(aes(color = species),
               width = 0.3, show.legend = FALSE) +
  geom_jitter(aes(color = species), alpha = 0.5, 
              show.legend = FALSE, 
              position = position_jitter(width = 0.2, seed = 0)) +
  scale_color_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "Epecie",
       y = "Longitud de la aleta en mm") +
  theme_bw() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20))


## ---- eval=FALSE--------------------------------------------
## ggplot(penguins) +
##   geom_point(mapping = aes(x = flipper_length_mm,
##                            y = body_mass_g,
##                            color = sex), size=3)+ theme_bw() +
##   theme(axis.text = element_text(size=20),
##         axis.title = element_text(size=20, face = "bold"),
##         legend.title = element_text(size=20)) +
##   guides(fill = guide_legend(title = "Sexo"))


## ---- out.height='60%', out.width='60%', fig.align='center', message=FALSE, warning=FALSE, echo=FALSE----
ggplot(penguins) +
  geom_point(mapping = aes(x=flipper_length_mm,y = body_mass_g,color = sex), size=3)+ theme_bw()+
  labs(x = "Longitud de la aleta en mm",
       y = "Masa corporal en gramos") +
  guides(fill = guide_legend(title = "Sexo")) +  
theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20, face = "bold"),
        legend.title = element_text(size=20)) 


## ---- message=FALSE, warning=FALSE--------------------------
penguins_long <- penguins %>% 
  tidyr::pivot_longer(
    cols = contains("mm"),
    names_to = "var", values_to = "val") %>% 
  tidyr::drop_na()


## ---- eval=FALSE--------------------------------------------
## penguins_long %>%
##   ggplot(aes(x = var, y = val, fill=var)) +
##   geom_boxplot() +
##   theme_bw() +
##   theme(legend.position="none") +
##   labs(x="", y="Medida en mm")


## ---- echo=FALSE, out.height='50%', out.width='50%', fig.align='center', message=FALSE, warning=FALSE----
penguins_long %>% 
  ggplot(aes(x = var, y = val, fill=var)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position="none") +
  labs(x="", y="Medida en mm") +
  theme(axis.text = element_text(size=20),
axis.title = element_text(size=20, face = "bold"),
legend.title = element_text(size=20)) 


## ---- out.height='40%', out.width='60%', fig.align='center'----
penguins %>%
  dplyr::count(species) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  ggplot() + geom_col(aes(x = prop, y = species)) +
  labs(x="Proporción", y="") + 
  theme(axis.text = element_text(size=20),
axis.title = element_text(size=20, face = "bold"),
legend.title = element_text(size=20)) 


## ---- echo=FALSE, out.width='80%',fig.align='center'--------
knitr::include_graphics("Imgs/extensiones.png")


## ---- eval=FALSE--------------------------------------------
## library(GGally)
## library(gapminder)
## gapminder %>% select(-country,-year) %>%
##   ggpairs(aes(color=continent))


## ---- echo=FALSE,message=FALSE, warning=FALSE, out.width='60%',fig.align='center'----
library(GGally)
library(gapminder)
gapminder %>% select(-country,-year) %>% 
  ggpairs(aes(color=continent))


## ---- echo=FALSE, out.width='50%', out.height='50%', fig.align='center'----
knitr::include_graphics("Imgs/violin.PNG")

